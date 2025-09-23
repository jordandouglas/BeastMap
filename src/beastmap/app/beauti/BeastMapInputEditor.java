package beastmap.app.beauti;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.ThreadedTreeLikelihood;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beastfx.app.inputeditor.BEASTObjectInputEditor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor.ExpandOption;
import beastfx.app.util.FXUtils;
import beastfx.app.inputeditor.ListInputEditor;
import beastmap.evolution.BranchMutationSampler;
import beastmap.evolution.PatternlessAlignment;
import beastmap.evolution.StochasticMapper;
import beastmap.logger.BranchSubstLogger;
import beastmap.logger.TypedTreeLogger;
import beastmap.logger.mut.SubstitutionSum;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import beast.base.core.Loggable;
import beastmap.logger.SampledSubstTreeLogger;
import beastmap.logger.SubstitutionSummer;



public class BeastMapInputEditor extends ListInputEditor {

	
	static HashMap<String, BranchMutationSampler> mapperMap = new HashMap<>();
	static HashMap<String, Logger> segmentedMap = new HashMap<>();
	static HashMap<String, Logger> substMap = new HashMap<>();
	
	public BeastMapInputEditor() {
		
	}
	
    public BeastMapInputEditor(BeautiDoc doc) {
        super(doc);
        Log.warning("BeastMapInputEditor 1");
    }
    
    
    @Override
	public Class<?> type() {
		return List.class;
	}

	@Override
	public Class<?> baseType() {
		return BranchMutationSampler.class;
	}

//	@Override
//	public Class<?> baseType() {
//		return Tree.class;
//	}
	
	
    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
    	
    	Log.warning("BeastMapInputEditor 2");
    	
    	super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    	
		
		pane = FXUtils.newVBox();
		
	    
	    // General instructions
		HBox label1 = FXUtils.newHBox();
		label1.getChildren().add(new Label("Select which loggers to include in the analysis."));
		HBox label2 = FXUtils.newHBox();
		label2.getChildren().add(new Label("\tSegmented tree loggers produce one branch segment each time the state changes (ideal for discrete geography or subsequences)."));
		HBox label3 = FXUtils.newHBox();
		label3.getChildren().add(new Label("\tBranch subst loggers count the number of events on each branch."));
		pane.getChildren().add(label1);
		pane.getChildren().add(label2);
		pane.getChildren().add(label3);
		
		
		// Prepare table
		TableView<LoggerSelection> table = new TableView<>();
		table.setEditable(true);
		
		// Prepare checkboxes
		CheckBox selecteAllCheckBox = new CheckBox();
	    	selecteAllCheckBox.setOnAction(
		        event -> {
		          event.consume();
		          table.getItems().forEach(
		        		  item -> item.setSelected(selecteAllCheckBox.isSelected())
		          );
		        });
	    	
    	
	    
	    // Table headers
	    TableColumn<LoggerSelection, String> nameCol = new TableColumn<>("Tree likelihood");
	    nameCol.setPrefWidth(250);
	    nameCol.setCellValueFactory(data -> data.getValue().nameProperty());
	    table.getColumns().add(nameCol);
	    TableColumn<LoggerSelection, String> loggerCol = new TableColumn<>("Logger type");
	    loggerCol.setPrefWidth(250);
	    loggerCol.setCellValueFactory(data -> data.getValue().loggerProperty());
	    table.getColumns().add(loggerCol);
	    
	    
	    // Check boxes
	    TableColumn<LoggerSelection, Boolean> selectedCol = new TableColumn<>("");
	    selectedCol.setGraphic(selecteAllCheckBox);
	    selectedCol.setSortable(false);
	    selectedCol.setPrefWidth(50);

	    selectedCol.setCellFactory( tc -> {
	    	CheckBoxTableCell cell =  new CheckBoxTableCell<>();
	    	cell.setEditable(true);
	    	List o = cell.getChildrenUnmodifiable();
	    	return cell;
	    });
	    selectedCol.setCellValueFactory( f -> f.getValue().selectedProperty());
	    selectedCol.setOnEditCommit(e->{
			Boolean newValue = e.getNewValue();
			LoggerSelection item = (LoggerSelection) e.getSource();
			item.setSelected(newValue);
		});
	    selectedCol.setEditable(true);
	    table.getColumns().add(selectedCol);


	    
	    
	    // Add table to page
	    pane.getChildren().add(table);
	    
  
		    
		
		// Get all tree likelihoods
		MCMC mcmc = (MCMC) doc.mcmc.get();
		CompoundDistribution posterior = (CompoundDistribution) mcmc.posteriorInput.get();
		List<GenericTreeLikelihood> likelihoods = new ArrayList<>();
		for (Distribution dist : posterior.pDistributions.get()) {
			if (dist.getID().equals("likelihood")) {
				CompoundDistribution cd = (CompoundDistribution) dist;
				for (Distribution likelihood : cd.pDistributions.get()) {
					if (likelihood instanceof GenericTreeLikelihood) {
						likelihoods.add((GenericTreeLikelihood)likelihood);
					}
				}
			}
		}
		
		
		// Trace log
		Logger traceLog = null;
		for (Logger logger : mcmc.loggersInput.get()) {
			if (logger.getID().equals("tracelog")) {
				traceLog = logger;
			}
		}
				
		
		// For each likelihood
		for (GenericTreeLikelihood likelihood : likelihoods) {

			
			
			String id = BeautiDoc.parsePartition(likelihood.getID());
			
			// Segmented logger
			BranchMutationSampler sampler;
			Logger segmentedTreeLogger, substTreeLogger;
			if (segmentedMap.containsKey(likelihood.getID())) {
				segmentedTreeLogger = segmentedMap.get(likelihood.getID());
				sampler = mapperMap.get(likelihood.getID());
				Log.warning("Found logger in map!");
			}else {
			
				// Get site model etc.
				SiteModelInterface siteModel = likelihood.siteModelInput.get();
				BranchRateModel.Base branchRateModel = likelihood.branchRateModelInput.get();
				TreeInterface tree = likelihood.treeInput.get();
				Alignment data = likelihood.dataInput.get();
				
				
				Log.warning(likelihood.getID() + " initialising data " + data);
				
				// Dataset
				PatternlessAlignment dataPatternless = new PatternlessAlignment();
				dataPatternless.initByName("data", data); // User data type??;
				
				
				Log.warning(likelihood.getID() + " initialising sampler");
				
				// Sampler
				sampler = new BranchMutationSampler();
				sampler.initByName("siteModel", siteModel, "branchRateModel", branchRateModel, "tree", tree, "burnin", 100000, "data", dataPatternless);
				
				Log.warning(likelihood.getID() + " initialising logger");
				
				
				// Loggable
				List<Loggable> segmentedLoggables = new ArrayList<>();
				TypedTreeLogger typedTreeLogger = new TypedTreeLogger();
				typedTreeLogger.initByName("sampler", sampler, "tree", tree);
				segmentedLoggables.add(typedTreeLogger);
				
				// Loggers
				segmentedTreeLogger = new Logger();
				String outfile = "beastmap.segmented." + id + ".trees";
				segmentedTreeLogger.initByName("fileName", outfile, "logEvery", 10000, "mode", "tree", "log", segmentedLoggables);
				
				// Add to MCMC chain
				//mcmc.loggersInput.get().add(segmentedTreeLogger);
				segmentedMap.put(likelihood.getID(), segmentedTreeLogger);
				mapperMap.put(likelihood.getID(), sampler);
				Log.warning("Adding logger to map!");
				
			}
			
			
			// Subst logger
			if (substMap.containsKey(likelihood.getID())) {
				substTreeLogger = substMap.get(likelihood.getID());
				Log.warning("Found logger in map!");
			}else {

				
				Tree tree = sampler.getTree();
				
				// Samplers
				List<BranchSubstLogger> substCounters = new ArrayList<>();
				SubstitutionSum sum = new SubstitutionSum();
				sum.initByName("sampler", sampler);
				substCounters.add(sum);
				
				
				
				// Loggable
				List<Loggable> substLoggables = new ArrayList<>();
				SampledSubstTreeLogger typedTreeLogger = new SampledSubstTreeLogger();
				typedTreeLogger.initByName("tree", tree, "sampler", substCounters);
				substLoggables.add(typedTreeLogger);
				
				// Loggers
				substTreeLogger = new Logger();
				String outfile = "beastmap." + id + ".trees";
				substTreeLogger.initByName("fileName", outfile, "logEvery", 10000, "mode", "tree", "log", substLoggables);
				
				// Add to MCMC chain
				//mcmc.loggersInput.get().add(segmentedTreeLogger);
				substMap.put(likelihood.getID(), substTreeLogger);
				
				
				
				// Sum them all - this removes posterior, likelihood, and prior from the tracelog
//				SubstitutionSummer sum_summer = new SubstitutionSummer();
//				sum_summer.initByName("counter", sum);
//				traceLog.loggersInput.get().add(sum_summer);
				
				
				Log.warning("Adding logger to map!");
				
			}
			
			
			// Add to table
			table.getItems().add(new LoggerSelection(likelihood.getID(), "Segmented tree logger", segmentedTreeLogger, mcmc));
			table.getItems().add(new LoggerSelection(likelihood.getID(), "Substitution count logger", substTreeLogger, mcmc));

			
		}
		
		

		
		getChildren().add(pane);
    	
    }
    
    
    // Include a logger in MCMC?
    public class LoggerSelection {
    	
    	private StringProperty name = new SimpleStringProperty();
    	private StringProperty loggerType = new SimpleStringProperty();
    	private MCMC mcmc;
		private Logger logger;
	    private BooleanProperty selected = new SimpleBooleanProperty();
	    
	    
	    public void setSelected(Boolean selected) { 
	    	this.selected.set(selected); 
			
	    	if (selected) {
	    		if (!this.mcmc.loggersInput.get().contains(logger)) {
	    			mcmc.loggersInput.get().add(logger);
	    		}
	    	}else {
	    		if (this.mcmc.loggersInput.get().contains(logger)) {
	    			mcmc.loggersInput.get().remove(logger);
	    		}
	    	}
	    	
	    }
	    public Boolean isSelected() { return selected.get(); }
	    public BooleanProperty selectedProperty() { return selected; }
	   
	    public void setName(String name) { this.name.set(name); }
	    public String getName() { return name.get(); }
	    public StringProperty nameProperty() { return name; }
	    
	    public void setLoggerType(String type) { this.loggerType.set(type); }
	    public StringProperty loggerProperty() { return loggerType; }

	    public LoggerSelection(String name, String loggerType, Logger logger, MCMC mcmc) {
	    	
	    	this.setName(name);
	    	this.setLoggerType(loggerType);
	    	this.logger = logger;
	    	this.mcmc = mcmc;
	    	
			setSelected(this.mcmc.loggersInput.get().contains(logger));
			selected.addListener((obs,  wasSelected,  isSelected) -> {
				setSelected(isSelected);
			});	
	    }
	    
	  }
    
    
    public void refreshPanel() {
       
    }
    
    public static boolean customConnector(BeautiDoc doc) {
    	
    	
    	Log.warning("BEASTMAP HERE");
    	return true;
    	
    }
    
}


