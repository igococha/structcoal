package beast.phylodynamics.model;

import org.jblas.DoubleMatrix;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;

public class RatesInputModel extends CalculationNode {
    public Input<RealParameter> migrationFromOutsideInput = new Input<RealParameter>(
    		"migrationFromOutside",
    		"offset of seasonality",
    		Input.Validate.REQUIRED);
    public Input<RealParameter> transmissionRateInput = new Input<RealParameter>(
    		"transmissionRate",
            "the rates of transmission between and within demes, the dimension needs to be n x n.",
            Validate.REQUIRED);
    public Input<RealParameter> migrationRateInput = new Input<RealParameter>(
    		"migrationRate",
            "the rates of migration between and within demes, the dimension needs to be n x n.",
            Validate.REQUIRED);    
    public Input<RealParameter> recoveryRateInput = new Input<RealParameter>(
    		"recoveryRate",
            "the number of all hosts at time of origin.",
            Input.Validate.REQUIRED);
    public Input<RealParameter> immunityLossInput = new Input<RealParameter>(
    		"immunityLoss",
            "the number of all hosts at time of origin.",
            Input.Validate.REQUIRED);    
    public Input<RealParameter> dimensionInput = new Input<RealParameter>(
    		"dimension",
            "number of demes/number of different populations (defaults to 2).",
            Input.Validate.REQUIRED);
    
    /*
     * Population Input Parameters
     */
    public Input<RealParameter> susceptibleFractionInput = new Input<RealParameter>(
    		"susceptible",
            "the number of susceptible hosts at time of origin (defaults to 0.5)");
    public Input<RealParameter> infectedFractionInput = new Input<RealParameter>(
    		"infected",
            "the number of infected hosts at time of origin (defaults to .1)");
    public Input<RealParameter> popSizeInput = new Input<RealParameter>(
    		"popSize",
            "the number of all hosts at time of origin (defaults to 3000).",
            Input.Validate.REQUIRED); 
    public Input<RealParameter> popScaleInput = new Input<RealParameter>(
    		"popScale",
            "always takes the n-times value for population size of the first population"); 
    public Input<BooleanParameter> migrationAsTransmissionInput = new Input<>(
    		"migrationAsTransmission",
            "weather to model state changes as transmission or migration"); 
    /*
     * Definition of accessible fields
     */
    
    protected DoubleMatrix migrationFromOutside;
    protected DoubleMatrix migrationRate;
    protected DoubleMatrix transmissionRate;
    protected DoubleMatrix recoveryRate;
    protected DoubleMatrix immunityLoss;
    protected DoubleMatrix popSize;
    protected DoubleMatrix sFraction;
    protected DoubleMatrix iFraction;
    protected int states;
    protected boolean symMig = false;
    protected boolean minTransRate = false;
    public boolean migrationAsTransmission = false;
    
    private boolean dirty = true;
    
	@Override
	public void initAndValidate() throws Exception {
		if(migrationAsTransmissionInput.get()!=null) 
			migrationAsTransmission = migrationAsTransmissionInput.get().getValue();
		updateMatrices();
	}
	
	public void update(){
		updateMatrices();
	}
	
	public void updateMatrices(){		
		states = (int) dimensionInput.get().getValue().doubleValue();
	    migrationFromOutside = DoubleMatrix.zeros(states);
	    migrationRate = DoubleMatrix.zeros(states,states);
	    transmissionRate = DoubleMatrix.zeros(states,states);
	    recoveryRate = DoubleMatrix.zeros(states);
	    immunityLoss = DoubleMatrix.zeros(states);
	    popSize = DoubleMatrix.zeros(states);
	    sFraction = DoubleMatrix.zeros(states);
	    iFraction = DoubleMatrix.zeros(states);
	    
	    updateMigrationTransmission();
	    
	    
		/*
		 * Migration from Outside (dampens internal oscillations)
		 */
		if(migrationFromOutsideInput.get() == null){
			for (int i = 0; i < states; i ++){
				migrationFromOutside.put(i, 0);
			}
		}else{
			RealParameter outsideMigration = migrationFromOutsideInput.get();
			if (outsideMigration.getDimension()==1){
				for (int i = 0; i < states; i ++)
					migrationFromOutside.put(i, migrationFromOutsideInput.get().getValue(i));
			}else if(outsideMigration.getDimension()==states){
				for (int i = 0; i < states; i ++)
					migrationFromOutside.put(i, migrationFromOutsideInput.get().getValue());
			}else {
				throw new IllegalArgumentException("The parameter migrationFromOutside has the wrong dimensionality");
			}
		}
		
    	/*
    	 * Recovery Rate Input. Can have either dim = states or dim = 1
    	 */
    	
    	RealParameter recRate = recoveryRateInput.get();
    	if (recRate.getDimension() == states){
    		for (int s = 0; s < states; s++)
    			recoveryRate.put(s,  recRate.getValue(s));
    	} else if (recRate.getDimension() == 1){
    		for (int s = 0; s < states; s++)
    			recoveryRate.put(s,  recRate.getValue());
    	} else {
    		throw new IllegalArgumentException("wrong number of elements in the recovery rate input"); 
    	}
    	
    	/*
    	 * Immunity Loss Input (wanning immunity). Can have either dim = states or dim = 1
    	 */
    	
    	RealParameter imLoss = immunityLossInput.get();
    	if (imLoss.getDimension() == states){
    		for (int s = 0; s < states; s++)
    			immunityLoss.put(s,  imLoss.getValue(s));
    	} else if (imLoss.getDimension() == 1){
    		for (int s = 0; s < states; s++)
    			immunityLoss.put(s,  imLoss.getValue());
    	} else {
    		throw new IllegalArgumentException("wrong number of elements in the immunity loss input"); 
    	}
    	
    	/*
    	 * Susceptible and Infected Fraction Input
    	 */
    	RealParameter S0 = susceptibleFractionInput.get();
    	RealParameter I0 = infectedFractionInput.get();
    	if(susceptibleFractionInput != null){
    		if (S0.getDimension() == 1)
    			for (int s = 0; s < states; s++) sFraction.put(s, S0.getValue());
    		else if (S0.getDimension() == states)
    			for (int s = 0; s < states; s++) sFraction.put(s, S0.getValue(s));
    		else
    			throw new IllegalArgumentException("wrong number of elements in the susceptible fraction input"); 
    	}else{
    		for (int s = 0; s < states; s++) sFraction.put(s, 0.5);
    	}
    	
    	if(infectedFractionInput != null){
    		if (I0.getDimension() == 1)
    			for (int s = 0; s < states; s++) iFraction.put(s, I0.getValue());
    		else if (S0.getDimension() == states)
    			for (int s = 0; s < states; s++) iFraction.put(s, I0.getValue(s));
    		else
    			throw new IllegalArgumentException("wrong number of elements in the infected fraction input"); 
    	}else{
    		for (int s = 0; s < states; s++) sFraction.put(s, 0.01);
    	}
    	
    	/*
    	 * Overall Population Size Input
    	 */
    	RealParameter pSize = popSizeInput.get();
    	if (popScaleInput.get()==null){
	    	if(pSize.getDimension() == states){
	    		for (int s = 0; s < states; s++)
	    			popSize.put(s, pSize.getValue(s));
	    			
	    	} else if (pSize.getDimension() == 1){
	    		for (int s = 0; s < states; s++)
	    			popSize.put(s, pSize.getValue());    		
	    	} else {
	    		throw new IllegalArgumentException("wrong number of elements in the population size input"); 
	    	}
	    }else{
	    	if(pSize.getDimension()!=1){
	    		throw new IllegalArgumentException("wrong number of elements in the population size input" + 
	    				"only one input value expected together with the population size scaler");
	    	}else{
	    		for (int s = 0; s < states; s++)
	    			popSize.put(sFraction,  pSize.getValue() * popScaleInput.get().getValue(s));
	    	}
	    }
    	
	}
	
	private void updateMigrationTransmission(){		
		/* 
		 * Transmission rate matrix Input. Can have two different dimensions, 
		 * 1.: dim = states -> transmission only allowed within deme or
		 * 2.: states * states -> transmission allowed within as well 
		 * as in between demes
		 */
    	int count = 0;
    	RealParameter tRate = transmissionRateInput.get();
    	/*
    	 * Transmission can only happen in the diagonal
    	 */
    	if (tRate.getDimension() == states){
    		minTransRate = true;	// Flag to indicate that transmission only happens in the diagonals
    		transmissionRate = DoubleMatrix.zeros(states,states);
    		for (int s = 0; s < states; s++)
    			transmissionRate.put(s,s, tRate.getArrayValue(s));
    	}
    	/*
    	 * Transmission can happen within as well as between demes
    	 */
    	else if(tRate.getDimension() == states * states){
    		for (int source=0; source < states; source++)
    			for (int sink=0; sink < states; sink++){
    				transmissionRate.put(source, sink, tRate.getArrayValue(count)); count++;
    			}
    	}
    	/*
    	 * All transmission rates are the same and transmission only
    	 * happens within deme
    	 */
    	else if(tRate.getDimension() == 1){
    		for (int s=0; s < states; s++)
    			transmissionRate.put(s, s, tRate.getArrayValue());
    	} else {
    		throw new IllegalArgumentException("wrong number of elements in the transmission rate input"); 
    	}		
		/*
		 * Migration Matrix Input, can either be of dimension states*(states-1) -> unsymetric
		 * or of dimension states * (states-1)/2 -> symetric rate matrix
		 */    	
		RealParameter migRateIn = migrationRateInput.get();
		count = 0;
		/*
		 * If migration is modeled as transmission, but the transmission
		 * rate Input is taken from the migration rate matrix in order
		 * to estimate only of diagonal elements
		 */
		if(migrationAsTransmission){
			/*
			 * If migration is modeled as transmission events go threw the loops that are meant
			 * to be for the non diagonal tranmission matrices. Also set the G matrix, respectively the
			 * matrix for the migration to 0
			 */
			minTransRate = false;
	    	migrationRate = DoubleMatrix.zeros(states,states);

			/*
			 * Rates from i to j are the same as rates from j to i
			 */
	    	if (migRateIn.getDimension() == (states * (states-1))/2){
	    		symMig = true;
	    		for (int source = 0; source < states; source++ ){
	    			for (int sink = source + 1; sink < states; sink++){
	    				transmissionRate.put(source, sink, migRateIn.getValue(count));
	    				transmissionRate.put(sink, source, migRateIn.getValue(count));
	    				count++;
	    			}
	    		}
	    	}
	    	/*
	    	 * All rates can be different
	    	 */
	    	else if (migRateIn.getDimension() == (states * (states-1))){
	    		for (int source = 0; source < states; source++ ){
	    			for (int sink = 0; sink < states; sink++){
	    				if(source!=sink){
	    					transmissionRate.put(source, sink, migRateIn.getValue(count));
	    					count++;
	    				}
	    			}
	    		}
	    	} 
	    	else {
	    		throw new IllegalArgumentException("wrong number of elements in the migration rate input"); 
	    	}
	    }
		/*
		 * Model migration events as events that are independent of the sink
		 * demes susceptible population, whereas when modelling threw transmission,
		 * the migration rate directly depends on the susceptible fraction Si/Ni
		 */
		else{
	    	if (migRateIn.getDimension() == (states * (states-1))/2){
//	    		symMig = true;
	    		for (int source = 0; source < states; source++ ){
	    			for (int sink = source + 1; sink < states; sink++){
	    				migrationRate.put(source, sink, migRateIn.getValue(count));
	    				migrationRate.put(sink, source, migRateIn.getValue(count));
	    				count++;
	    			}
	    		}
	    	}
	    	else if (migRateIn.getDimension() == (states * (states-1))){
	    		for (int source = 0; source < states; source++ ){
	    			for (int sink = 0; sink < states; sink++){
	    				if(source!=sink){
	    					migrationRate.put(source, sink, migRateIn.getValue(count));
	    					count++;
	    				}
	    			}
	    			migrationRate.put(source, source, 0.0); // fill the diagonal
	    		}
	    	} 
	    	else {
	    		throw new IllegalArgumentException("wrong number of elements in the migration rate input"); 
	    	}
	    }
	}
    
    public int getStates(){
    	return states;
    }
    
    public DoubleMatrix getMigrationFromOutside(){
    	return migrationFromOutside;
    }
    
    public DoubleMatrix getMigrationRate(){
    	return migrationRate;
    }
    
    public DoubleMatrix getTransmissionRate(){
    	return transmissionRate;
    }
    
    public DoubleMatrix getRecoveryRate(){
    	return recoveryRate;
    }
    
    public DoubleMatrix getImmunityLoss(){
    	return immunityLoss;
    }
    
    public DoubleMatrix getPopSize(){
    	return popSize;
    }
    
    public DoubleMatrix getSusceptibleFraction(){
    	DoubleMatrix S = new DoubleMatrix(states);
    	for (int s = 0; s < states; s++)
    		S.put(s,  sFraction.get(s) * popSize.get(s));
    	return S;
    }
    
    public DoubleMatrix getInfectedFraction(){
    	DoubleMatrix I = new DoubleMatrix(states);
    	for (int s = 0; s < states; s++)
    		I.put(s,  iFraction.get(s) * popSize.get(s));
    	return I;
    }
    
    public DoubleMatrix getSFraction(){
    	return sFraction;
    }
    
    public DoubleMatrix getIFraction(){
    	return iFraction;
    }    
    
    public boolean getSymetricMigration(){
    	return symMig;
    }
    
    public boolean getMinimalTransmission(){
    	return minTransRate;
    }
    
    /*
     *  Check if any of the Parameter Input Values has changed
     *  If any of those values has changed, return true to show the
     *  EpiModel, that recalculation is required
     */   
    public boolean isDirty(){
    	if (transmissionRateInput.isDirty()) return true;
    	if (migrationRateInput.isDirty()) return true;   
    	if (recoveryRateInput.isDirty()) return true;
    	if (immunityLossInput.isDirty()) return true;
    	if (susceptibleFractionInput.isDirty()) return true;
    	if (infectedFractionInput.isDirty()) return true;
    	if (popSizeInput.isDirty()) return true;  
    	if (popScaleInput.isDirty()) return true;
    	return false;
    }
}
