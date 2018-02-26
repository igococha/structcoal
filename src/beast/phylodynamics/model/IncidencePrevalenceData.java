package beast.phylodynamics.model;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.datatype.DataType.Base;

public class IncidencePrevalenceData extends Base {
    public Input<String> timeSeriesInput = new Input<String>("timeSeries",
    		"mapping of time series of one deme to dates", Validate.REQUIRED);

    private List<DoubleMatrix> prevalence = new ArrayList<DoubleMatrix>();
    private List<Double> timePoint = new ArrayList<Double>();
    
    @Override
	public void initAndValidate() throws Exception {
        String[] prev = timeSeriesInput.get().split(",");
        String[] split;
    	DoubleMatrix addPrev = new DoubleMatrix(prev[0].split("=").length-1);
        for (int i = 0; i < prev.length; i++){
        	split = prev[i].split("=");
        	addPrev = new DoubleMatrix(split.length-1);
        	timePoint.add(Double.parseDouble(split[0]));        	
        	for (int j = 1; j < split.length; j++){
        		addPrev.put(j-1, Double.parseDouble(split[j]));
        	}
        	prevalence.add(addPrev);
        }
    }
    
	
	public double getTimeSeries(int state) {
        return timePoint.get(state);
    }
	
    public DoubleMatrix getNumberInfected(int state) {
        return prevalence.get(state);
    }
    
    public int getSeriesLength() {
        return prevalence.size();
    }

    public String getTypeDescription() {
        return "incidence Prevalence Data";
    }


}
