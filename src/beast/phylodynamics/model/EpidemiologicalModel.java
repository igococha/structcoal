package beast.phylodynamics.model;

import java.util.ArrayList;
import java.util.Collections;

import org.jblas.DoubleMatrix;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;


/**
 * @author Nicola Felix Mueller
 */
@Description("Deterministic compartmental epidemiological model with structure" +
             "Tracks deterministic pop trajectory and F, G and Y matrices")
public class EpidemiologicalModel extends CalculationNode {

	/*
	 * Input specifying the rates Object, which contains all the informations
	 * about the initial and estimated rates, such as migration, transmission
	 * or recovery rates
	 */
	public Input<RatesInputModel> ratesInput = new Input<RatesInputModel>("ratesInput",
			"Input of migration and transmission Rates",Validate.REQUIRED);	
    public Input<RealParameter> timeStepInput = new Input<RealParameter>("timeStep",
            "timeStep (defaults to 2).",Input.Validate.REQUIRED);    

    /*
     * Input specifying the origin, meaning the point at which to start with
     * simulating the S(E)IR trajectories as well as the f- g- and y-series 
     */
    public Input<RealParameter> originInput = new Input<RealParameter>("origin",
            "timepoint of origin.",Input.Validate.REQUIRED);

    /*
     * Inputs that are specific to the seasonal forcing model.   
     */
    public Input<RealParameter> seasonalStrengthInput = new Input<RealParameter>(
    		"seasonalStrength",
    		"how strong the seasonal patterns are, higher means less seasonal");  
    public Input<RealParameter> seasonalOffsetInput = new Input<RealParameter>(
    		"seasonalOffset",
    		"offset of the season");
    public Input<RealParameter> wavelengthInput = new Input<RealParameter>(
    		"wavelength",
    		"the wavelength of the oscillations");
    
    /*
     * Inputs specific to the exponential growth model, which is only the parameter
     * defining when a epidemic was first introduced into a deme
     */
    public Input<RealParameter> intitalIntroductionInput = new Input<RealParameter>(
    		"intitalIntroduction",
    		"when the first individualk entered the population");
    
    /*
     * Boolean Inputs used to define which model to use. Would technically not 
     * be necessary but is nice to have
     */    
    public Input<BooleanParameter> isSeasonalInput = new Input<BooleanParameter>(
    		"isSeasonal",
    		"define weather to use a seasonal forcing model");
    public Input<BooleanParameter> isExponentialInput = new Input<BooleanParameter>(
    		"isExponential",
    		"define weather to use a exponential growth model");
    public Input<BooleanParameter> isConstantInput = new Input<BooleanParameter>(
    		"isConstant",
    		"define weather to use a constant population model");

	
    
    /*
     * Parameters that are needed independent of the population model used
     */
    public boolean dirty;
    public boolean reject = false;
    public boolean diagTrans = false;    
    
    /*
     * TimeSeries Object that stores all the F, G, Y and S I R ArraLists
     */
    public TimeSeries timeS;
    
    /*
     * ArrayList storing all the times between which euler integration is
     * performed
     */
    public ArrayList<Double> integrationTimes;
   
    public int states;    
    protected DoubleMatrix migrationFromOutside;
    protected DoubleMatrix transmissionRate;
    public DoubleMatrix migrationRate;
    protected DoubleMatrix recRate;
    protected DoubleMatrix imLoss;
    
    /*
     * Parameters specific to seasonal forcing
     */
    private DoubleMatrix seasonalStrength;
    private DoubleMatrix seasonalOffset;
    
    /*
     * Parameters specific to the exponential growth model
     */
    private DoubleMatrix intitalInfected;
 
   
    /*
     * Boolean that store the population model used for 
     * calculating the F G & Y matrices over time. At least
     * one of those has to be true 
     */
    public boolean seasonal = false;
    public boolean exponential = false;
    public boolean constant = false;
    public boolean movingFraction = false;
    
    private double origin;
    private double timeStep;
    private double wavelength = 365;
    	
	@Override
	public void initAndValidate() throws Exception {
		diagTrans = ratesInput.get().minTransRate;
		states = ratesInput.get().getStates();
		timeStep = timeStepInput.get().getValue();
		
		if(originInput.get() != null) origin  = originInput.get().getValue();
		
		/*
		 * Get the population model used as Input
		 */
		if (isSeasonalInput.get() != null) seasonal = isSeasonalInput.get().getValue();
		if (isConstantInput.get() != null) constant = isConstantInput.get().getValue();
		if (isExponentialInput.get() != null) exponential = isExponentialInput.get().getValue();		
		if (wavelengthInput.get() != null) wavelength = wavelengthInput.get().getValue();
		
		timeS = new TimeSeries();
	}
	
	public boolean update(){
		
		timeS.setBack(states);		
    	boolean reject = false;
        
        final double endTime = 0.0;

        // Update the transmission and Migration Rate Matrices
        RatesInputModel rates = ratesInput.get();
        rates.update();         
        
        transmissionRate = rates.getTransmissionRate(); 
        migrationRate = rates.getMigrationRate();
        migrationFromOutside = rates.getMigrationFromOutside();
        recRate = rates.getRecoveryRate();
        imLoss = rates.getImmunityLoss();
        
        timeS.setS(rates.getSusceptibleFraction());
        timeS.setI(rates.getInfectedFraction());
        timeS.setN(rates.getPopSize());
        timeS.setR(timeS.getN().sub(timeS.getS().add(timeS.getI())));
        timeS.setPS(timeS.getPS());
        timeS.setPI(timeS.getPI());
        		
        integrationTimes = new ArrayList<Double>();    	
		
		double time = origin;
		
		while (time > endTime){
			integrationTimes.add(time);
			time -= timeStep;
		}		
		integrationTimes.add(endTime);
	
		timeS.addS(timeS.getS());
		timeS.addI(timeS.getI());
		timeS.addR(timeS.getN().sub(timeS.getS().add(timeS.getI())));
			
		/*
		 * Calculate using a seasonal model
		 * The seasonal Parameters need to be specified 
		 * as well as the boolean seasonal must be true
		 * otherwise an Error will be returned
		 */		
		if (seasonal){
			int t = 0;
			int dur = integrationTimes.size();
			updateSeasonalParam();
			
			do{
				reject = deltaSeasonal(t, dur, timeStep);
				t++;
			}while(integrationTimes.get(t-1)>0 && !reject);		

		}
		/*
		 * Calculate the trajectories using an exponential
		 * growth model. The parameter initialInfected
		 * must be specified as well as the boolean
		 * exponential needs to be set to true in order
		 * to work properly
		 */
		if (exponential){
			int t = 0;
			int dur = integrationTimes.size();
			boolean[] start = new boolean[states];
			updateInitialTimePoint();
			
			do{
				/*
				 * Add an initial individual at an estimated time point
				 * which is stored in the DoubleMatrix initialInfected 
				 */
				for (int s = 0; s < states; s++){
					if (!start[s] && intitalInfected.get(s) < integrationTimes.get(s)){
						timeS.setI(s, 1.0);
						start[s] = true;
					}
				}
				reject = deltaExponential(t, dur, timeStep);
				t++;
			}while(integrationTimes.get(t-1)>0 && !reject);			
		}
		/*
		 * fill the arrays with all the same values if
		 * a constant population mode is chosen
		 */
		if (constant){
			
			/*
			 * Initialize F G & Y Matrices to not have to get
			 * them in every loop
			 */
			DoubleMatrix F = rates.getTransmissionRate().mulColumnVector(timeS.getI()).mulRowVector(timeS.getS().div(timeS.getN()));
			DoubleMatrix G = rates.getMigrationRate().mulColumnVector(timeS.getI());
			DoubleMatrix Y = timeS.getI();

			/*
			 * Fill the timeSeries list with the data, as everything is supposed to be
			 * Independent of time, only the first entry is used
			 */
			timeS.addF(F);
			timeS.addG(G);
			timeS.addY(Y);
			timeS.addS(timeS.getS());
			timeS.addI(timeS.getI());
			timeS.addR(timeS.getR());
		}

		/*
		 * Reverse the array Lists in order to go from
		 * forward in time to backwards in time
		 */
		timeS.reverse();
		Collections.reverse(integrationTimes);

		return reject;		
	}
	
	/*
	 * Calculation of seasonal Trajectories as well as the updating of the seasonal Parameters
	 */
	protected boolean deltaSeasonal(int t, int dtTimeCount, double timeStep){
		
		DoubleMatrix currentSeason = new DoubleMatrix(states);
		DoubleMatrix dS = new DoubleMatrix(states);
		DoubleMatrix dI = new DoubleMatrix(states);
		DoubleMatrix dR = new DoubleMatrix(states);
		
		for (int d = 0; d < states; d++){
			currentSeason.put(d,Math.sin(2*Math.PI/wavelength*(integrationTimes.get(t)-seasonalOffset.get(d))));	// Value of sinus at time t 
		}
	
//		imLoss.print();
		DoubleMatrix season = seasonalStrength.add(currentSeason).div(seasonalStrength.add(1));			// add seasonality term
		DoubleMatrix seasRates = new DoubleMatrix();
		seasRates.copy(transmissionRate);
		for (int d = 0; d < states; d++){
			seasRates.put(d, d, transmissionRate.get(d,d)*season.get(d)); 								// multiply row one of transmission rates with element 1 of seasonality vector etc
		}
		
//		seasRates = transmissionRate.mulRowVector(season);
		
		/*
		 *  Calculate the change in the individual fractions (susceptibles, infected, recovered)
		 */		
		DoubleMatrix Transmission = seasRates.mulColumnVector(timeS.getI()).mulRowVector(timeS.getS().div(timeS.getN())); 
		DoubleMatrix Recovery = timeS.getI().mul(recRate);		
		DoubleMatrix ImmunityLoss = timeS.getR().mul(imLoss);
		DoubleMatrix Migration = migrationRate.mulColumnVector(timeS.getI()); 
			
		
		/*
		 *  Update dS/dt dI/dt and dR/dt
		 */
		dS.copy(ImmunityLoss.sub(Transmission.columnSums()));
		dS = dS.add(Migration.rowSums());
		dS = dS.sub(Migration.columnSums());
		
		dI.copy(Transmission.columnSums().sub(Recovery));
		dI = dI.add(Migration.columnSums());
		dI = dI.sub(Migration.rowSums());
	
		dR.copy(Recovery.sub(ImmunityLoss)); 
		

		
		// Update F, G, and Y
		timeS.addF(Transmission);
		timeS.addG(Migration);            
	    timeS.addY(timeS.getI());
	    if (t < (dtTimeCount - 1)) {
	    	// Update state variables
	    	final double dt = integrationTimes.get(t) - integrationTimes.get(t+1);
	    	
	    	timeS.setS(timeS.getS().add(dS.mmul(dt)));
	    	timeS.setI(timeS.getI().add(dI.mmul(dt)));
	    	timeS.setR(timeS.getR().add(dR.mmul(dt)));
	    	timeS.setN(timeS.getS().add(timeS.getI().add(timeS.getR())));
	    		    			
	    	
	    	if (timeS.getS().min() < 0 || timeS.getI().min() < 0){
	    		return true;
	    	}
           	timeS.addS(timeS.getS());
           	timeS.addI(timeS.getI());
           	timeS.addR(timeS.getR());
   	
	    }

	    return false;
	}

	private void updateSeasonalParam(){

		seasonalStrength = new DoubleMatrix(states);
	    seasonalOffset = new DoubleMatrix(states);

		for(int s = 0; s <states; s++)
			seasonalStrength.put(s, seasonalStrengthInput.get().getValue(s));
		RealParameter seasOffset = seasonalOffsetInput.get();
		if(seasOffset.getDimension() == states){
			for(int s = 0; s <states; s++)
				seasonalOffset.put(s, seasOffset.getValue(s));
		} else{
			for(int s = 0; s <states; s++)
				seasonalOffset.put(s, seasOffset.getValue());					
		}
	}
	
	/*
	 * Calculation of the exponential trajectories as well as updating of
	 * the initial time point of the first individual
	 */
	protected boolean deltaExponential(int t, int dtTimeCount, double TimeStep){
		
		DoubleMatrix Transmission = transmissionRate.mulColumnVector(timeS.getI()).mulRowVector(timeS.getS().div(timeS.getN()));
		DoubleMatrix Recovery = timeS.getI().mul(recRate);		
		DoubleMatrix ImmunityLoss = timeS.getR().mul(imLoss);
		DoubleMatrix Migration = migrationRate.mulColumnVector(timeS.getS());
		
		DoubleMatrix dS = new DoubleMatrix(states);
		DoubleMatrix dI = new DoubleMatrix(states);
		DoubleMatrix dR = new DoubleMatrix(states);
		/*
		 *  Update dS/dt dI/dt and dR/dt
		 */
		dS.copy(ImmunityLoss.sub(Transmission.columnSums()));
		dI.copy(Transmission.columnSums().sub(Recovery));
		dR.copy(Recovery.sub(ImmunityLoss));
		
		timeS.addF(Transmission);
		timeS.addG(Migration);
		timeS.addY(timeS.getI());
		
		if (t < (dtTimeCount - 1)) {
		    
	    	// Update state variables
	    	final double dt = integrationTimes.get(t) - integrationTimes.get(t+1);
	    		   	
	    	timeS.setS(timeS.getS().add(dS.mmul(dt)));
	    	timeS.setI(timeS.getI().add(dI.mmul(dt)));
	    	timeS.setR(timeS.getR().add(dR.mmul(dt)));
	    	timeS.setN(timeS.getS().add(timeS.getI().add(timeS.getR())));
	    		    			
	    	if (timeS.getPS().min() < 0|| timeS.getPI().min() < 0){
	    		return true;
	    	}
           	timeS.addS(timeS.getS());
           	timeS.addI(timeS.getI());
           	timeS.addR(timeS.getR());
	    }
		
		return false;
	}
	
	private void updateInitialTimePoint(){
		intitalInfected = new DoubleMatrix(states);
		for (int s = 0; s < states; s++)
			intitalInfected.put(s, intitalIntroductionInput.get().getValue(s));
	}
	
	/*
	 * something else
	 */
	protected boolean delta(int i, int t, double s){
		return false;
	}
	
	public int getTrajLength(){
    	return timeS.size();
    }

    public DoubleMatrix getNumS(int t) {
    	return timeS.getS(t);
    }
    
	public DoubleMatrix getNumI(int t) {
		return timeS.getI(t);
    }
    
	public DoubleMatrix getNumR(int t) {
		return timeS.getR(t);
    }
	
	public double getNumSDeme(int t, int d) {
		return timeS.getS(t).get(d);
    }
    
	public double getNumIDeme(int t, int d) {
		return timeS.getI(t).get(d);
    }
    
	public double getNumRDeme(int t, int d) {
		return timeS.getR(t).get(d);
    }
	
	public double getTime(int t) {
		return integrationTimes.get(t);
    }
    
    /**
     * @param t time interval to return
     * @param k state to return Y for
    */
    public double getY(int t, int k) {
    	if (constant)
    		return timeS.getY(0).get(k);
    	else
    		return timeS.getY(t).get(k);
    }
    
    public DoubleMatrix getY(int t) {
    	if (constant)
    		return timeS.getY(0);
    	else
    		return timeS.getY(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in F
     * @param l column in F
    */
    public double getF(int t, int k, int l) {
    	if (constant)
    		return timeS.getF(0).get(k,l);
    	else
    		return timeS.getF(t).get(k,l);
    }
    
    public DoubleMatrix getF(int t) {
    	if (constant)
    		return timeS.getF(0);
    	else
    		return timeS.getF(t);
    }
    
    /**
     * @param t time interval to return
     * @param k row in G
     * @param l column in G
    */
    public double getG(int t, int k, int l) {
    	if (constant)
    		return timeS.getG(0).get(k,l);
    	else
    		return timeS.getG(t).get(k,l);
    }
    
    public DoubleMatrix getG(int t) {
	    if (constant)
			return timeS.getG(0);
		else
			return timeS.getG(t);
    }
    
    /**
     * 
     * @param t
     * @param k
     * @param l
     * @return
     */
    public double getM(int t, int k, int l) {
    	if (constant)
    		return timeS.getM(0).get(k,l);
    	else
    		return timeS.getM(t).get(k,l);
    }
    
    public DoubleMatrix getM(int t) {
	    if (constant)
			return timeS.getM(0);
		else
			return timeS.getM(t);
    }
    
    public Boolean movingFraction(){
    	return movingFraction;
    }
    
    /**
     * @return integrationTimes in reverse order (backwards in time)
     */
    public ArrayList<Double> getIntegrationTimes() {
    	return integrationTimes;
    }
    
    public double getOrigin(){
    	return origin;
    }
    
    public void setOrigin(double origin){
    	this.origin = origin;
    }
    
    /*
     * CalculationNode interface
     */

    @Override
    public boolean requiresRecalculation() {
        dirty = true;
        return true;
    }

    @Override
    public void restore() {
        dirty = true;
        super.restore();
    }
    
    public boolean isDirty() {
    	if (seasonalStrengthInput.isDirty()) return true; 
    	if (seasonalOffsetInput.isDirty()) return true;
    	if (intitalIntroductionInput.isDirty()) return true;
    	return ratesInput.get().isDirty();
    }
	
	
}
