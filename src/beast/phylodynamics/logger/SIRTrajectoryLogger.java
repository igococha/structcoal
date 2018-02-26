package beast.phylodynamics.logger;

import java.io.PrintStream;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.phylodynamics.model.EpidemiologicalModel;


/**
 * @author Nicola Mueller
 */
@Description("Logs the migration rates with source and sink")
public class SIRTrajectoryLogger extends BEASTObject implements Loggable {
    public Input<EpidemiologicalModel> epiModelInput = new Input<EpidemiologicalModel>(
    		"epiModel",
            "the rates of migration between and within demes, the dimension needs to be n x n.",
            Validate.REQUIRED);
    
    protected int states;
    
    @Override
    public void initAndValidate() throws Exception {

    }
	

    @Override
    public void init(PrintStream out) throws Exception {
    	epiModelInput.get().update();
    	for (int i = 0; i < epiModelInput.get().getTrajLength(); i++){
			if ( i < epiModelInput.get().getTrajLength()-1){
				out.print("t_" + i + "\t");
			}else{
				out.print("t_" + i);
			}
    	}
    	out.print("\n");
    	out.print("0" + "\t");
    	out.print("0" + "\t");
    	for (int i = 1; i < epiModelInput.get().getTrajLength(); i++){
				if ( i < epiModelInput.get().getTrajLength()-1){
					out.print(epiModelInput.get().getTime(i) + "\t");
				}else{
					out.print(epiModelInput.get().getTime(i));
				}
	    	}
    }
    
    @Override
    public void log(int nSample, PrintStream out) {
    	epiModelInput.get().update();

        //Susceptible Fraction;
    	for (int d = 0; d < epiModelInput.get().states; d++){
    		if (d==0){
    		}
    		else{
    			out.print(nSample + "\t");
    		}
    		for (int i = 0; i < epiModelInput.get().getTrajLength(); i++){
				if ( i < epiModelInput.get().getTrajLength()-1){
					out.format("%.2f\t",epiModelInput.get().getNumSDeme(i, d));
				}else{
					out.format("%.2f",epiModelInput.get().getNumSDeme(i, d));
				}
	    	}
	    	out.print("\n");
    	}
    	
    	//Infected Fraction;
    	for (int d = 0; d < epiModelInput.get().states; d++){
   			out.print(nSample + "\t");

    		for (int i = 0; i < epiModelInput.get().getTrajLength(); i++){
				if ( i < epiModelInput.get().getTrajLength()-1){
					out.format("%.2f\t",epiModelInput.get().getNumIDeme(i, d));
				}else{
					out.format("%.2f",epiModelInput.get().getNumIDeme(i, d));
				}				
	    	}
	    	out.print("\n");
    	}
    	
    	//Recovered Fraction;
    	for (int d = 0; d < epiModelInput.get().states; d++){
   			out.print(nSample + "\t");
    		for (int i = 0; i < epiModelInput.get().getTrajLength(); i++){
				if ( i < epiModelInput.get().getTrajLength()-1){
					out.format("%.2f\t",epiModelInput.get().getNumRDeme(i, d));
				}else{
					out.format("%.2f",epiModelInput.get().getNumRDeme(i, d));
				}		
	    	}
    		if(d < (epiModelInput.get().states-1)){
    			out.print("\n");
    		}
    	}  	
    }
    
    @Override
    public void close(PrintStream out) {
    }

}
