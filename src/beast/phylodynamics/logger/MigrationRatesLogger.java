package beast.phylodynamics.logger;

import java.io.PrintStream;

import org.jblas.DoubleMatrix;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.phylodynamics.model.RatesInputModel;


/**
 * @author Nicola Mueller
 */
@Description("Logs the migration rates with source and sink")

public class MigrationRatesLogger extends BEASTObject implements Loggable {
	public Input<RatesInputModel> ratesInput = new Input<RatesInputModel>(
			"ratesInput",
			"Input of migration and transmission Rates",
			Validate.REQUIRED);	
	public Input<BooleanParameter> migrationInput = new Input<BooleanParameter>(
			"migration",
			"Input of migration");
    
    protected int states;
    private boolean trans = true;
    
    @Override
    public void initAndValidate() throws Exception {
    	states = ratesInput.get().getStates();
    	if(migrationInput.get()!=null)
    		trans = migrationInput.get().getValue();
    }
	
	/*
     * Loggable interface
     */
    @Override
    public void init(PrintStream out) throws Exception {
    	if(ratesInput.get().migrationAsTransmission){
    		for (int source = 0; source < states; source++)
    			for (int sink = 0; sink < states; sink++)
	    				out.print("transmission" + source + "to" + sink + "\t");
    	}else{
	    	if (ratesInput.get().getSymetricMigration()){
	    		for (int source = 0; source < states; source++)
	    			for (int sink = source+1; sink < states; sink++)
	    				if(source != sink) 				
		    				out.print("symMig_" + source + "to" + sink + "\t");
	    	}else{
	    		for (int source = 0; source < states; source++)
	    			for (int sink = 0; sink < states; sink++)
	    				if(source != sink)
		    				out.print("mig_" + source + "to" + sink + "\t");
	    	}
	    }
    }
    
    @Override
    public void log(int nSample, PrintStream out) {
        //out.format("%g\t", dt);
    	if (ratesInput.get().migrationAsTransmission){
	    	DoubleMatrix transRates = ratesInput.get().getTransmissionRate();
    		for (int source = 0; source < states; source++){
    			for (int sink = 0; sink < states; sink++){
    					out.print(transRates.get(sink,source) + "\t");
    			}
    		}	    	
    	}else{
	    	DoubleMatrix migRates = ratesInput.get().getMigrationRate();
	    	if (ratesInput.get().getSymetricMigration())
	    		for (int source = 0; source < states; source++)
	    			for (int sink = source; sink < states; sink++)
	    				if(sink!=source)
	    					out.print(migRates.get(sink,source) + "\t");
	    	
	    	/*
	    	 * Print migration Rates Backward in time
	    	 */
	    	if(!ratesInput.get().getSymetricMigration())
	    		for (int source = 0; source < states; source++)
	    			for (int sink = 0; sink < states; sink++)
	    				if(sink!=source)
	    					out.print(migRates.get(sink,source) + "\t");
	   	}
    }
    
    @Override
    public void close(PrintStream out) {
    }

}
