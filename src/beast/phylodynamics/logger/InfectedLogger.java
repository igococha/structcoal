package beast.phylodynamics.logger;

import java.io.PrintStream;

import org.jblas.DoubleMatrix;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.phylodynamics.model.RatesInputModel;


/**
 * @author Nicola Mueller
 */
@Description("Logs the migration rates with source and sink")

public class InfectedLogger extends BEASTObject implements Loggable {
	public Input<RatesInputModel> ratesInput = new Input<RatesInputModel>(
			"ratesInput",
			"Input of migration and transmission Rates",
			Validate.REQUIRED);	
    
    protected int states;
    
    @Override
    public void initAndValidate() throws Exception {
    	states = ratesInput.get().getStates();    	
    }
	
	/*
     * Loggable interface
     */
    @Override
    public void init(PrintStream out) throws Exception {
    	for (int s = 0; s < states; s++)
    		out.print("I" + s + "\t");
    }
    
    @Override
    public void log(int nSample, PrintStream out) {
        //out.format("%g\t", dt);
    	DoubleMatrix inf = ratesInput.get().getInfectedFraction();
    	for (int s = 0; s < states; s++)
    		out.print(inf.get(s) + "\t");
    }
    
    @Override
    public void close(PrintStream out) {
    }

}
