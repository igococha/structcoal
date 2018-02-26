package beast.phylodynamics.logger;

import java.io.PrintStream;

import org.jblas.DoubleMatrix;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.phylodynamics.epidemiology.StructCoalVolzDensity;

public class ProbabilisticTreeLogger extends BEASTObject implements Loggable {
    public Input<StructCoalVolzDensity> coalescentDensityInput = new Input<StructCoalVolzDensity>(
    		"coalescentDensity",
    		"tree over which to calculate a prior or likelihood");
    public Input<RealParameter> dimensionInput = new Input<RealParameter>(
    		"dimension",
            "number of demes/number of different populations (defaults to 2).",
            Input.Validate.REQUIRED);


    
    protected int states;
    protected StructCoalVolzDensity strCoal;
    
    @Override
    public void initAndValidate() throws Exception {
    	double tmp = dimensionInput.get().getValue();
    	states = (int) (tmp);

    }
	
	
	/*
     * Loggable interface
     */
  @Override
  public void init(PrintStream out) throws Exception {
	  	for (int i = 0; i < (coalescentDensityInput.get().nrSamples-1); i++){
	  		for (int j = 0; j < states; j++){
	  			int nodeNr = coalescentDensityInput.get().nrSamples + i + 1; 
	  			out.print("n_" + nodeNr +"d_" + j + "\t");
	  		}
	  	}
  }

  @Override
  public void log(int nSample, PrintStream out) {
//	try {
//		coalescentDensityInput.get().calculateLogP();
//	} catch (Exception e) {
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	}

	  DoubleMatrix[] stateProbs = coalescentDensityInput.get().getStateProbabilities();
	  	for (int i = 0; i < stateProbs.length; i++){
	  		for (int j = 0; j < states; j++){
	  			out.format("%.3f\t",stateProbs[i].get(j));
	  			}
	  	}
  }

  @Override
  public void close(PrintStream out) {
  }
  
}
