package beast.phylodynamics.logger;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;

import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.StateNode;
import beast.core.Input.Validate;
import beast.core.parameter.Parameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.phylodynamics.epidemiology.StructCoalVolzDensity;


/**
 * @author Nicola Felix Mueller
 */
public class StructuredTreeLogger extends Tree implements Loggable {
	public Input<StructCoalVolzDensity> volzDensityInput = new Input<StructCoalVolzDensity>(
			"volzDensity",
			"A deterministic epidemiological model",
			Input.Validate.REQUIRED);
    public Input<Tree> treeInput = new Input<Tree>("tree", "tree to be logged", Validate.REQUIRED);

    public Input<BranchRateModel.Base> clockModelInput = new Input<BranchRateModel.Base>("branchratemodel", "rate to be logged with branches of the tree");
    public Input<List<Function>> parameterInput = new Input<List<Function>>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());
    public Input<Boolean> maxStateInput = new Input<Boolean>("maxState", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Boolean> substitutionsInput = new Input<Boolean>("substitutions", "report branch lengths as substitutions (branch length times clock rate for the branch)", false);
    public Input<Integer> decimalPlacesInput = new Input<Integer>("dp", "the number of decimal places to use writing branch lengths and rates, use -1 for full precision (default = full precision)", -1);
    public Input<TraitSet> typeTraitInput = new Input<>(
            "typeTrait", "Type trait set.");
    
    boolean someMetaDataNeedsLogging;
    boolean substitutions = false;
    boolean takeMax = true;

    private DecimalFormat df;
    private String type;
    private boolean traitInput = false;
	
	
    @Override
    public void initAndValidate() throws Exception {
    	
        if (typeTraitInput.get() != null) traitInput = true;

    	
        if (parameterInput.get().size() == 0 && clockModelInput.get() == null) {
        	someMetaDataNeedsLogging = false;
        	return;
            //throw new Exception("At least one of the metadata and branchratemodel inputs must be defined");
        }
    	someMetaDataNeedsLogging = true;
    	// without substitution model, reporting substitutions == reporting branch lengths 
        if (clockModelInput.get() != null) {
        	substitutions = substitutionsInput.get();
        }
       
        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();
        	
        }

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }
        
        type = volzDensityInput.get().getType();
    }

    @Override
    public void init(PrintStream out) throws Exception {
    	treeInput.get().init(out);
    }

    @Override
    public void log(int nSample, PrintStream out) {
        // make sure we get the current version of the inputs
        Tree tree = (Tree) treeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i < metadata.size(); i++) {
        	if (metadata.get(i) instanceof StateNode) {
        		metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        	}
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
//        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadata, branchRateModel));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }

    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel) {
        if (maxStateInput.get() != null){
        	takeMax = maxStateInput.get();
        	
        }
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
        if (!node.isLeaf()) {
        	if (!takeMax){	        
		        buf.append("[&" + type + "prob={");
		        DoubleMatrix stateProbs = volzDensityInput.get().getStateProb(node.getNr());
		        
		        for (int i = 0 ; i < volzDensityInput.get().states-1; i++){
		        	buf.append(String.format("%.3f", stateProbs.get(i)));
		        	buf.append(",");
		        }
		        buf.append(String.format("%.3f", stateProbs.get(volzDensityInput.get().states-1)));
		        buf.append("}");
	        
		        buf.append(",max" + type + "=");
		        buf.append(String.format("%d", stateProbs.argmax() ));
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");
		        DoubleMatrix stateProbs = volzDensityInput.get().getStateProb(node.getNr());

		        buf.append(String.format("%d", stateProbs.argmax() ));
		        buf.append(']');        		
        	}
        }else{
			String sampleID = node.getID();
			String[] splits = sampleID.split("_");
			int sampleState;
			
			if(traitInput){				
				sampleState = (int) typeTraitInput.get().getValue(node.getID());
			}
			else{
				sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
			}
			if (!takeMax){
    	        
		        buf.append("[&" + type + "prob={");
	
		        for (int i = 0 ; i < volzDensityInput.get().states-1; i++){
		        	if (sampleState != i) buf.append(String.format("0.0"));
		        	if (sampleState == i) buf.append(String.format("1.0"));
	            	buf.append(",");
		        }

	        	if (sampleState != volzDensityInput.get().states-1) buf.append(String.format("0.0"));
	        	if (sampleState == volzDensityInput.get().states-1) buf.append(String.format("1.0"));
		        
		        buf.append("}");
		        buf.append(",max" + type + "=");

		        buf.append(String.format("%d", sampleState ));
//		        buf.append(']'); 
		        buf.append(']');
        	}else{
		        buf.append("[&max" + type + "=");

		        buf.append(String.format("%d", sampleState ));
		        buf.append(']');        		
        	}
        }
        
        buf.append(":");
        if (substitutions) {
            appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, node.getLength());
        }
        return buf.toString();
    }


    @Override
    public void close(PrintStream out) {
    	treeInput.get().close(out);
    }

	
}
