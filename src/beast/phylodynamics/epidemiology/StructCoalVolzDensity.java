package beast.phylodynamics.epidemiology;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.coalescent.IntervalType;
import beast.phylodynamics.model.EpidemiologicalModel;


/**
 * @author David Rasmussen
 */

@Description("Calculates the probability of a beast.tree using under an unstructured SIR model using the framework of Volz (2009).")
public class StructCoalVolzDensity extends StructuredTreeDistribution implements Loggable {

	public Input<EpidemiologicalModel> epiModelInput = new Input<EpidemiologicalModel>(
			"epiModel",
			"A deterministic epidemiological model",
			Input.Validate.REQUIRED);
	
    public Input<TraitSet> typeTraitInput = new Input<>(
            "typeTrait", "Type trait set.");
    
    public Input<BooleanParameter> useEulerInput = new Input<>(
    		"useEuler",
            "wether to use euler integration to update the linestateProbabilities(default = false)"); 
    
    public Input<IntegerParameter> EulerStepsInput = new Input<>(
    		"EulerSteps",
            "in how many steps euler integration is to perform (Default 1)"); 
	
	public int samples;
	public int nrSamples;
	public DoubleMatrix[] stateProbabilities;
	public DoubleMatrix[] stateProbabilities1;
	public DoubleMatrix[] stateProbabilities2;
    
    public int states;
    
    private boolean needsUpdate = true;
    private boolean traitInput = false;
    private boolean useEuler = false;
    private int EulerSteps = 1;
    
    // Set up for lineage state probabilities
    ArrayList<Integer> activeLineages;
    ArrayList<DoubleMatrix> lineStateProbs;
    
    @Override
    public void initAndValidate() throws Exception {
    	treeIntervalsInput.get().calculateIntervals();
 	
        if (treeIntervalsInput.get() == null)
            throw new Exception("Expected treeIntervals to be specified");
        
        stateProbabilities = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];        
        stateProbabilities1 = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];        
        stateProbabilities2 = new DoubleMatrix[treeIntervalsInput.get().getSampleCount()];        
        nrSamples = treeIntervalsInput.get().getSampleCount() + 1;

        states = epiModelInput.get().states;
        
        if (typeTraitInput.get() != null) traitInput = true;
        if (useEulerInput.get() != null) useEuler = useEulerInput.get().getValue();
        if (EulerStepsInput.get() != null ) EulerSteps = (int) (EulerStepsInput.get().getValue());
        
        calculateLogP();
    }

    public double calculateLogP() throws Exception {
    	boolean reject = epiModelInput.get().update();
    	treeIntervalsInput.get().calculateIntervals();
    	
    	// if origin is being sampled it must be constrained to be older than the root of the tree.
    	if (epiModelInput.get().getOrigin() < treeIntervalsInput.get().getTotalDuration()) {
    		if (epiModelInput.get().constant){
    			System.out.println("");
    			System.out.println("origin is set to the tree hight*1.01");
    			System.out.println("");
    			epiModelInput.get().setOrigin(treeIntervalsInput.get().getTotalDuration()*1.01);
    			epiModelInput.get().update();
    		}else{
	    		System.out.println("origin is too low");
	            logP = Double.NEGATIVE_INFINITY;
	            System.out.println("Coal logP is NEG INepiModelInput.get().update()F");
	            return logP;
    		}
        }

    	// if maximal count exceeded or number too small
        if (reject) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }
        
        // Set up for lineage state probabilities
        activeLineages = new ArrayList<Integer>(); 
        lineStateProbs = new ArrayList<DoubleMatrix>();
        
        // Compute likelihood at each integration time and tree event starting at final sampling time and moving backwards
        ArrayList<Double> integrationTimes = epiModelInput.get().getIntegrationTimes();							//in reverse time
        int currTreeInterval = 0; 																				// what tree interval are we in?
        double currTime = integrationTimes.get(0); 																// current time (since present)
        double nextIntegrationTime; 																			// time of next integration step
        double nextIntervalTime = treeIntervalsInput.get().getInterval(0); 										// time next tree interval starts
        final int intervalCount = treeIntervalsInput.get().getIntervalCount();
        
        logP = 0;
        int t = 0;
        
        do {        	
        	nextIntegrationTime = integrationTimes.get(t+1);													// Update Integration time
        	while (nextIntervalTime <= nextIntegrationTime) {													// while there are still events to take place before the next integration step
        		/*
        		 * compute contribution of last interval to likelihood
        		 */
	        	final double duration = nextIntervalTime - currTime;
	        	if (duration > 0) {
	        		double add = -computeLambdaSum(t) * duration;
	        		if (add > 0)
	        		{
	        			System.out.println("lambda sum is larger than one:\t" + add + "\t" + duration);
	        			epiModelInput.get().getF(t).print();
	        			epiModelInput.get().getY(t).print();
	        		}	
	        		logP += add;
        			/*
        			 * newly added here to update line states until a coalescent event happens
        			 * if not here the migration rates might be biased in regions were a lot of
        			 * coalescent events happen
        			 */
	  				boolean negInf = updateLineProbs(t, duration); 
	  				if (negInf)
	  					logP = Double.NEGATIVE_INFINITY;
	        	}                
	               
	        	/*
	        	 *  compute contribution of event to likelihood
	        	 */
	        	if (treeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.COALESCENT) {	//CHECK IF INDEX IS CORRECT HERE
	        		final double yTotal = computeYTotal(t);
	        		if (yTotal < 1.0) {
	        			logP = Double.NEGATIVE_INFINITY;
	        			return logP;
	        		}
	        		
	        		final double pairCoalRate = computeLambdaPair(t, currTreeInterval);
	        		if (pairCoalRate >= 0) {
	        			logP += Math.log(pairCoalRate);
	        		}
	        		updateLineProbsCoalEvent(t, currTreeInterval);												// Set parent lineage state probs and remove children
		       	}
	       		
	        	/*
	        	 * add new lineage
	        	 */
	       		if (treeIntervalsInput.get().getIntervalType(currTreeInterval) == IntervalType.SAMPLE) {
	       			addLineages(currTreeInterval);
	       		}
	       		
	       		currTime += duration;
	       		currTreeInterval++; 
	       		
	       		if (currTreeInterval < intervalCount) {
	       			nextIntervalTime = currTime + treeIntervalsInput.get().getInterval(currTreeInterval);
        		} else {
	       			nextIntervalTime = Double.POSITIVE_INFINITY;
	       			currTreeInterval--; //stay in last interval
	       		}

        	}   
        	
	        final double duration = nextIntegrationTime - currTime;
	    	if (duration > 0) {    			
	   			// compute contribution of last time step to likelihood
  				logP += -computeLambdaSum(t) * duration; // if A is greater than Y???
  				boolean negInf = updateLineProbs(t, duration); 
  				if (negInf)
  					logP = Double.NEGATIVE_INFINITY;
    		}		    	
    		currTime += duration;
    		t++;
    		
        }while(integrationTimes.get(t-1) <= treeIntervalsInput.get().getTotalDuration());
	        
        //Compute likelihood of remaining tree intervals (coal events occuring before origin)
        if (Double.isInfinite(logP))logP = Double.NEGATIVE_INFINITY;
        return logP;
   	
    }
    
    /*
     * Should be correct
     */
    private double computeLambdaPair(int t, int currTreeInterval) throws Exception {
		
		List<Node> coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
    	if (coalLines.size() > 2) {
			throw new Exception("Unsupported coalescent at non-binary node");
		}
    	if (coalLines.size() < 2) {
    		System.out.println();
    		System.out.println("WARNING: Less than two lineages found at coalescent event!");
    		System.out.println();
    		return Double.NaN;
		}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			return 0.0;
		}
		double lambda = 0.0;
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
//		System.out.println(t);
//		epiModelInput.get().getY(t).print();
        for (int k = 0; k < states; k++) { 
        	double Yk = epiModelInput.get().getY(t,k);
            if ( Yk > 0.0) {
            	for (int l = 0; l < states; l++){
                	double Yl = epiModelInput.get().getY(t,l);
                	if (Yl > 0.0){
		            	final double popCoalRate = epiModelInput.get().getF(t,k,l) / (Yk*Yl);
						final double pairCoalRate = popCoalRate * 
							(lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(l)
								+ lineStateProbs.get(daughterIndex1).get(l) * lineStateProbs.get(daughterIndex2).get(k));
						
						if (!Double.isNaN(pairCoalRate)) lambda += pairCoalRate;
						
                	}
            	}
            }
        }
		
    	return lambda;
    }
    
    /*
     * Should be correct, matrix stuff not checked
     */
    private double computeLambdaSum(int t) {
    	
		int lineCount = activeLineages.size();
		double lambdaSum = 0.0; 										// holds the sum of the pairwise coalescent rates over all lineage pairs
		
		if (lineCount < 2) return lambdaSum;											// Zero probability of two lineages coalescing	


		/*
		 * Sum over line pairs (scales better with numbers of demes and lineages)
		 */	
		DoubleMatrix popCoalRates;
		/*
		 * If the transmission matrix has only elements in the diagonal
		 * do fancy matrix calculation epiModelInput.get().diagTrans
		 */
    	if(epiModelInput.get().diagTrans){
			popCoalRates = new DoubleMatrix(states);
			for (int k = 0; k < states; k++){
				final double Y = epiModelInput.get().getY(t,k);
				if(Y>0)
					popCoalRates.put(k, (epiModelInput.get().getF(t,k,k) / (Y*Y)));
			}
    		
			DoubleMatrix sumStates = DoubleMatrix.zeros(states);
			DoubleMatrix diagElements = DoubleMatrix.zeros(states);
			for (int linI = 0; linI < lineCount; linI++) {
				sumStates = sumStates.add(lineStateProbs.get(linI));
				diagElements = diagElements.add(lineStateProbs.get(linI).mul(lineStateProbs.get(linI)));
			}
			DoubleMatrix M = sumStates.mul(sumStates);
			lambdaSum = M.sub(diagElements).mul(popCoalRates).sum();
			
		/*
		 * if the transmission matrix has off diagonal elements, use the standard
		 * way of calculating the interval contribution	
		 */
    	}else{    		
			popCoalRates = new DoubleMatrix(states, states);
			for (int k = 0; k < states; k++) {
				for (int l = 0; l < states; l++) {
					final double Yk = epiModelInput.get().getY(t,k);
					final double Yl = epiModelInput.get().getY(t,l);
					if (Yk > 0 & Yl > 0) {
						popCoalRates.put(k,l, (epiModelInput.get().getF(t,k,l) / (Yk*Yl)));
					} else {
						popCoalRates.put(k,l, 0);
					}
				}	
			}
			for (int linI = 0; linI < lineCount; linI++) {					
				for (int linJ = linI+1; linJ < lineCount; linJ++) {


					//Sum over line states
					for (int k = 0; k < states; k++) {
						for (int l = 0; l < states; l++) {
							if (k==l){	
								final double pairCoalRate = popCoalRates.get(k,l) * 
									(lineStateProbs.get(linI).get(k) * lineStateProbs.get(linJ).get(l)
										+ lineStateProbs.get(linI).get(l) * lineStateProbs.get(linJ).get(k));
								if (!Double.isNaN(pairCoalRate)) lambdaSum += pairCoalRate;
							}

						}
					}
				}	
			}
    	}
    	return lambdaSum;    	
    }
    
    private boolean updateLineProbs(int t, double dt) {    	
		/*
		 * Set up Q matrix for lineage state transition rates. The 
		 * G Matrix is filled with the migration rates times the source
		 * Population forward in time.
		 */
    	/*
    	 * If the transmission Rate flag minTRate is true, meaning
    	 * that transmissions can only happen within a deme, than
    	 * states changed are only calculated using the G matrix
    	 * & the Y matrix(Infected Population)
    	 */    	
//    	if(epiModelInput.get().diagTrans){
//    		Q = epiModelInput.get().getG(t).divColumnVector(epiModelInput.get().getY(t));  // Take only the migration Rate as Input & not the G-matrix as lineages move independently of population dynamics 
       	/*
    	 * Else there is contribution of the transmission of lineages in deme k
    	 * to deme l, so the contribution to state change by transmission needs 
    	 * to be taken into account
    	 */
//    	}else{
    		/*
    		 * Get the different measures: The A matrix is the summed probability 
    		 * over every lineage of beeing in a state i. The Y matrix is the 
    		 * number of infected in state i. The mul matrix is (Y - A) ./ Y 
    		 */
		DoubleMatrix A = getLineStateSum();
		DoubleMatrix Y = epiModelInput.get().getY(t);
		DoubleMatrix mul = Y.sub(A).div(Y);  
		DoubleMatrix M = new DoubleMatrix();
		
		/*
		 * If there is more state probability in a deme than individuals
		 * then return the false flag indicating an error
		 */
		if(mul.min()<0)
			return true;
		if (Double.isNaN(mul.get(0)))
			return true;
			
		/*
		 * Calculate the rate at which states change. State changes in time t
		 * happen from k to l. The total rate at which a state changes depends on
		 * the migration (G-matrix) and the transmission form deme k to l (F-matrix)
		 */
    	DoubleMatrix Q = DoubleMatrix.zeros(states,states);
   	
		for (int k = 0; k < states; k++) 
		{
			double rowSum = 0.0;
			for (int l = 0; l < states; l++) 
			{
				if (k != l) 
				{																        				// off-diagonal
					if (Y.get(k) > 0) 
					{
						final double rateOutByBirth = epiModelInput.get().getF(t,l,k) * mul.get(l) / Y.get(k);       // Birth in other deme (Backwards in time)
						final double rateOutByMigration = epiModelInput.get().getG(t,l,k) / (Y.get(k)); 								// State change of lineages (Backwards in time)
						final double totalRate = rateOutByBirth + rateOutByMigration;
						Q.put(k, l, totalRate);
						rowSum += totalRate;
					} 
					else
					{																// diagonal
						Q.put(k, l, 0.0);
					}
				}
			}
			Q.put(k,k,-rowSum);
		}
		/*
		 * Update lineage state probabilities for all active lineages. The change in state
		 * probability is calculated as: i) Q the state change matrix per time is multiplied 
		 * with time ii) Q is multiplied with the current state probabilities such that
		 * dP = Q*dt*s has the elements [q11 * s1 , q12 * s1...; q21 *s2...].
		 * The changes in state probabilities is than the sum of columns.  
		 */
		/** Test to see how much worse euler does then Matrix Exponentiation
		if(true)
		{	DoubleMatrix C = new DoubleMatrix();
			C.copy(M);
			DoubleMatrix newP = MatrixFunctions.expm(C);
			for (int lin = 0; lin < activeLineages.size(); lin++) {
				DoubleMatrix oldProbs = lineStateProbs.get(lin);			
				DoubleMatrix probs = vecMat(newP, oldProbs);			
				probs = probs.div(probs.sum());
				probs.print();
				DoubleMatrix dP = C.muliColumnVector(oldProbs); //returns change in state probs
				lineStateProbs.get(lin).add(dP.columnSums()).print(); //add gained/lost probability mass to state prob vector				
				System.out.println("");
			}		
		}
		*/
		if(!useEuler) // calculate linestate probs using euler integration
		{
			DoubleMatrix newP = matrixExponential(Q, dt);
			for (int lin = 0; lin < activeLineages.size(); lin++) {
				DoubleMatrix oldProbs = lineStateProbs.get(lin);			
				DoubleMatrix probs = vecMat(newP, oldProbs);			
				/*
				 * divide the state probabilities by the sum of state probabilities
				 * in order to avoid errors due to numerical errors as small as they
				 * might be
				 */
				probs = probs.div(probs.sum());
				lineStateProbs.set(lin, probs);
			}
		}else
		{
			DoubleMatrix C = new DoubleMatrix();		
			double dt_euler = dt/EulerSteps;
			C = Q.mul(dt_euler);
			for (int eul = 0; eul < EulerSteps; eul++){			
				for (int lin = 0; lin < activeLineages.size(); lin++) {
					DoubleMatrix oldProbs = lineStateProbs.get(lin);	
					DoubleMatrix dP = C.mulColumnVector(oldProbs); //returns change in state probs
					DoubleMatrix probs = oldProbs.add(dP.columnSums());		
					lineStateProbs.set(lin, probs.div(probs.sum())); //add gained/lost probability mass to state prob vector		
				}
			}
		}


		
		/*
		 * Everything seems to have worked correctly
		 */
		return false;
		
    }
    
    private void updateLineProbsCoalEvent(int t, int currTreeInterval) throws Exception {
    	List<Node> coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
    	if (coalLines.size() > 2) {
			throw new Exception("Unsupported coalescent at non-binary node");
		}
		int daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
		int daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
		
		if (daughterIndex1 == -1 || daughterIndex2 == -1){
	    	treeIntervalsInput.get().swap();
	    	coalLines = treeIntervalsInput.get().getLineagesRemoved(currTreeInterval);
	    	daughterIndex1 = activeLineages.indexOf(coalLines.get(0).getNr());
	    	daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr());
	    	if (daughterIndex1 == -1 || daughterIndex2 == -1)
				throw new Exception("Active lineages does not contain coalescing lineages");
		}
		List<Node> parentLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		if (parentLines.size() > 1) throw new Exception("Unsupported coalescent at non-binary node");
		
		
		//Add parent to activeLineage and initialize parent's state prob vector
		Node parentNode = parentLines.get(0);
		activeLineages.add(parentNode.getNr());
		DoubleMatrix pVec = new DoubleMatrix(states);

		//Compute parent lineage state probabilities in pVec
		DoubleMatrix coalRates = DoubleMatrix.zeros(states, states);
		double lambda;

		for (int k = 0; k < states; k++) {
            for (int l = 0; l < states; l++) {            	
            	final double Yk = epiModelInput.get().getY(t,k);
            	final double Yl = epiModelInput.get().getY(t,l);
            	if (Yk > 0 && Yl > 0)
            		lambda = epiModelInput.get().getF(t,k,l) / (Yk*Yl) 
            			* (lineStateProbs.get(daughterIndex1).get(k) * lineStateProbs.get(daughterIndex2).get(l)
            				+ lineStateProbs.get(daughterIndex1).get(l) * lineStateProbs.get(daughterIndex2).get(k));
            	else 
            		lambda = 0.0;
            	
            	coalRates.put(k, l, lambda);
            }
        }	
		/*
		 * pVec := new state probabilities of the lineage dt after coalescing		
		 */
		pVec = coalRates.rowSums().div(coalRates.sum());
//		boolean isnan = false;
//		for (int i = 1; i < states; i++)
//			if (Double.isNaN(pVec.get(i)))
//				isnan = true;
//		if (isnan){
//			System.out.println("NaN at node");
//			epiModelInput.get().getG(t).print();
//			epiModelInput.get().getF(t).print();
//			epiModelInput.get().getY(t).print();
//			pVec.print();
//			System.out.println(computeLambdaPair(t, currTreeInterval));
//			System.out.println(lineStateProbs.get(daughterIndex1).get(0));
//			System.out.println(lineStateProbs.get(daughterIndex1).get(1));
//			System.out.println(lineStateProbs.get(daughterIndex2).get(0));
//			System.out.println(lineStateProbs.get(daughterIndex2).get(1));
//			if (computeLambdaPair(t, currTreeInterval) > 0.0000001)
//				System.exit(1);
//			else 
//				System.out.println("");
//			
//		}
		
		stateProbabilities[parentNode.getNr() - nrSamples] = (pVec);
        stateProbabilities1[parentNode.getNr() - nrSamples] = lineStateProbs.get(daughterIndex1);
        stateProbabilities2[parentNode.getNr() - nrSamples] = lineStateProbs.get(daughterIndex2);
        lineStateProbs.add(pVec);
        
		//Remove daughter lineages
		activeLineages.remove(daughterIndex1);
		lineStateProbs.remove(daughterIndex1);
		
		daughterIndex2 = activeLineages.indexOf(coalLines.get(1).getNr()); //index may have changed after removal of first daughter
		activeLineages.remove(daughterIndex2);
		lineStateProbs.remove(daughterIndex2);
    }
    
    private void addLineages(int currTreeInterval) {
		List<Node> incomingLines = treeIntervalsInput.get().getLineagesAdded(currTreeInterval);
		if(traitInput){
			/*
			 * If there is a typeTrait given as Input the model will take this
			 * trait as states for the taxons
			 */		
			for (Node l : incomingLines) {				
				activeLineages.add(l.getNr());
				int sampleState = (int) typeTraitInput.get().getValue(l.getID());
				DoubleMatrix sVec = DoubleMatrix.zeros(states);
				sVec.put(sampleState, 1.0);			
				lineStateProbs.add(sVec);
			}			
		}else{		
			/*
			 * If there is no trait given as Input, the model will simply assume that
			 * the last value of the taxon name, the last value after a _, is an integer
			 * that gives the type of that taxon
			 */
			for (Node l : incomingLines) {
				activeLineages.add(l.getNr());
				String sampleID = l.getID();
				String[] splits = sampleID.split("_");
				int sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
				DoubleMatrix sVec = DoubleMatrix.zeros(states);
				sVec.put(sampleState, 1.0);			
				lineStateProbs.add(sVec);
			}	
		}
    }
    
    private DoubleMatrix getLineStateSum() {    	
    	DoubleMatrix sumAk = DoubleMatrix.zeros(states);
   		for (int lin = 0; lin < activeLineages.size(); lin++)
   			sumAk = sumAk.add(lineStateProbs.get(lin));
    	return sumAk;    	
    }
    
    private double computeYTotal(int t) {
    	return epiModelInput.get().getY(t).sum();
    }
    
    public DoubleMatrix getStateProb(int nr){
    	return stateProbabilities[nr - nrSamples];
    }
    public DoubleMatrix getStateProb1(int nr){
    	return stateProbabilities1[nr - nrSamples];
    }    
    public DoubleMatrix getStateProb2(int nr){
    	return stateProbabilities2[nr - nrSamples];
    }    
    public DoubleMatrix[] getStateProbabilities(){
    	return stateProbabilities;
    }
    
    public String getType(){
    	if (typeTraitInput.get()!=null) return typeTraitInput.get().getTraitName();
    	else return "type";
    }
    
    @Override
    protected boolean requiresRecalculation() {
    	return true;
    }
    
    /*
     * Quickly added vector times matrix 
     */
    private DoubleMatrix vecMat(DoubleMatrix mat, DoubleMatrix vec){
    	DoubleMatrix prod = DoubleMatrix.zeros(states);
    	for (int i = 0; i < states; i++)
    		for (int j = 0; j < states; j++)
    			prod.put(i, prod.get(i) + vec.get(j)*mat.get(j,i));
    	return prod;
    }
    
    private DoubleMatrix matrixExponential(DoubleMatrix Q, double dt){
    	return MatrixFunctions.expm(Q.mul(dt));
    }
}