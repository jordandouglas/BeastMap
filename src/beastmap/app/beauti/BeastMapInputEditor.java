package beastmap.app.beauti;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.atomic.AtomicBoolean;

import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.parser.PartitionContext;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import beastfx.app.inputeditor.ListInputEditor;
import beastmap.evolution.BranchMutationSampler;
import beastmap.evolution.PatternlessAlignment;
import beastmap.logger.BranchSubstLogger;
import beastmap.logger.TypedTreeLogger;
import beastmap.logger.mut.*;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.ReadOnlyObjectWrapper;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.scene.control.CheckBox;
import javafx.scene.control.Label;
import javafx.scene.control.TableCell;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.control.cell.CheckBoxTableCell;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;
import starbeast3.tree.SpeciesTree;
import beast.base.core.Loggable;
import beastmap.logger.SampledSubstTreeLogger;
import beastmap.logger.StochasticMapProperty;
import beastmap.logger.SubstitutionSummer;
import beastmap.logger.SubstitutionSummerPerBranch;






public class BeastMapInputEditor extends ListInputEditor {

	
	public static Class[] termsToLog = new Class[] { SubstitutionSum.class, SynonymousSubstSum.class, NonSynonymousSubstSum.class, 
														NucleotideTransitionCounter.class, NucleotideTransversionCounter.class, 
														AminoAcidClassChanges.class, AminoAcidClassRemains.class
										  }; // TotalSize.class NetSize.class, FromToSubstSum.class,
	
	
	
	
	private SimpleIntegerProperty burnin = new SimpleIntegerProperty(100000);
	
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
	
	
	
	protected static TableView<LoggerSelection> table;
	
    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption, boolean addButtons) {
    	
    	
    	super.init(input, beastObject, itemNr, isExpandOption, addButtons);
    	
    	
    	
    	
    	
    	MCMC mcmc = (MCMC) doc.mcmc.get();
    	List<BranchMutationSampler> mappersAll = getAllMappers();
    	Log.warning("Found " + mappersAll.size() + " stochastic mappers in mcmc");
		
		pane = FXUtils.newVBox();
		
	    
	    // General instructions
		HBox label1 = FXUtils.newHBox();
		label1.getChildren().add(new Label("Select which loggers to include in the analysis."));
		HBox label2 = FXUtils.newHBox();
		label2.getChildren().add(new Label("\t•Segmented tree loggers produce one branch segment each time the state changes (ideal for discrete geography). Warning: tree files may grow very large."));
		HBox label3 = FXUtils.newHBox();
		label3.getChildren().add(new Label("\t•Branch substitution loggers count the number of events on each branch, but do not record the times."));
		pane.getChildren().add(label1);
		pane.getChildren().add(label2);
		pane.getChildren().add(label3);
		
		
		pane.getChildren().add(prepareBurninTextbox(mappersAll));
		
		
	    // Add table to page
		table = this.prepareTable();
	    pane.getChildren().add(table);
	    
		    
		
		// Get all tree likelihoods
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
				
		
		
		// Remove any unused likelihoods from the cache
		for (BranchMutationSampler mapper : mappersAll) {
			
			String id = mapper.getID();
			id = id.replaceAll("BeastMap.StochasticMapper.", "");
			
			Log.warning(id);
			
			boolean likelihoodActive = false;
			
			for (PartitionContext partition : doc.partitionNames) {
				
				Log.warning(partition.partition);
				
				if (partition.partition.equals(id)) {
					likelihoodActive = true;
					break;
				}
			}

			if (!likelihoodActive) {
				Log.warning("Removing " + id);
				
				String segmentedTreeLoggerID =  "BeastMap.SegmentedTreeLogger." + id;
				String substTreeLoggerID =  "BeastMap.SubstTreeLogger." + id;
				String mapperID =  "BeastMap.StochasticMapper." + id;
				
				removeWithID(mcmc.loggersInput.get(), substTreeLoggerID);
				removeWithID(mcmc.loggersInput.get(), segmentedTreeLoggerID);
				
				doc.unregisterPlugin(doc.pluginmap.get(segmentedTreeLoggerID));
				doc.unregisterPlugin(doc.pluginmap.get(substTreeLoggerID));
				doc.unregisterPlugin(doc.pluginmap.get(mapperID));
				
			}
			
		}
		
		
	
		
		
		// For each likelihood
		for (GenericTreeLikelihood likelihood : likelihoods) {

			
			
			String id = BeautiDoc.parsePartition(likelihood.getID());
			String segmentedTreeLoggerID =  "BeastMap.SegmentedTreeLogger." + id;
			String substTreeLoggerID =  "BeastMap.SubstTreeLogger." + id;
			String mapperID =  "BeastMap.StochasticMapper." + id;
			
			
			
			// Mapper/sampler
			BranchMutationSampler sampler;
			if (inputContainsID(mappersAll, mapperID)) {
				sampler = (BranchMutationSampler) getInputWithID(mappersAll, mapperID);
				Log.warning("Found mapper in mcmc! " + mapperID);
			}else {
				
				
				Alignment data = likelihood.dataInput.get();
				
				
				// Dataset
				PatternlessAlignment dataPatternless = new PatternlessAlignment();
				dataPatternless.initByName("data", data); // User data type??;
				dataPatternless.setID("BeastMap.PatternlessAlignment." + id);
				
				
				// Sampler
				sampler = new BranchMutationSampler();
				sampler.initByName("likelihood", likelihood, "data", dataPatternless, "burnin", burnin.get());
				sampler.setID(mapperID);
				
				
				
				//doc.registerPlugin(dataPatternless);
				doc.registerPlugin(sampler);
				
				
				mappersAll.add(sampler);
				Log.warning("Making new mapper " + mapperID);
				
			}
			
			
			// Segmented logger
			Logger segmentedTreeLogger, substTreeLogger;
			if (inputContainsID(mcmc.loggersInput.get(), segmentedTreeLoggerID)) {
				segmentedTreeLogger = (Logger) getInputWithID(mcmc.loggersInput.get(), segmentedTreeLoggerID);
				Log.warning("Found logger in map! " + segmentedTreeLoggerID);
				
			}else {
			
			
			
				
				// Loggable
				List<Loggable> segmentedLoggables = new ArrayList<>();
				TypedTreeLogger typedTreeLogger = new TypedTreeLogger();
				typedTreeLogger.initByName("sampler", sampler, "likelihood", likelihood);
				typedTreeLogger.setID("BeastMap.TypedTreeLogger." + id);
				segmentedLoggables.add(typedTreeLogger);
				
				// Loggers
				segmentedTreeLogger = new Logger();
				String outfile = "beastmap.segmented." + id + ".trees";
				segmentedTreeLogger.initByName("fileName", outfile, "logEvery", 10000, "mode", "tree", "log", segmentedLoggables);
				segmentedTreeLogger.setID(segmentedTreeLoggerID);
				
				
				doc.registerPlugin(typedTreeLogger);
				doc.registerPlugin(segmentedTreeLogger);
				
				Log.warning("Making new mapper " + segmentedTreeLoggerID);

				
			}
			
			
			// Subst logger
			if (inputContainsID(mcmc.loggersInput.get(), substTreeLoggerID)) {
				substTreeLogger = (Logger) getInputWithID(mcmc.loggersInput.get(), substTreeLoggerID);
				Log.warning("Found logger in map! " + substTreeLoggerID);
			} else {


				
				substTreeLogger = getTreeLoggerForLikelihood((Tree) likelihood.treeInput.get());
				
				if (substTreeLogger == null) {
					
					// Loggable
					List<Loggable> substLoggables = new ArrayList<>();
					SampledSubstTreeLogger typedTreeLogger = new SampledSubstTreeLogger();
					typedTreeLogger.initByName("likelihood", likelihood); //, "sampler", substCounters);
					typedTreeLogger.setID("BeastMap.SampledSubstTreeLogger." + id);
					substLoggables.add(typedTreeLogger);
					
					// Loggers
					substTreeLogger = new Logger();
					String outfile = "beastmap." + id + ".trees";
					substTreeLogger.initByName("fileName", outfile, "logEvery", 10000, "mode", "tree", "log", substLoggables);
					substTreeLogger.setID(substTreeLoggerID);
					

					doc.registerPlugin(typedTreeLogger);
					doc.registerPlugin(substTreeLogger);
					
					
				}
				

				
				Log.warning("Making new mapper " + substTreeLoggerID);
				
			}
			
			
			
			SampledSubstTreeLogger substLogger = null;
			for (BEASTObject obj : substTreeLogger.loggersInput.get()) {
				if (obj instanceof SampledSubstTreeLogger) {
					substLogger = (SampledSubstTreeLogger) obj;
				}
			}
			
			
			List<SubstLoggerSelection> loggers = getSubstSumLoggers(sampler, id, substLogger, likelihood.dataInput.get().getDataType());
			
			
			
			// Add to table
			table.getItems().add(new LoggerSelection(id, "Segmented tree logger", segmentedTreeLogger, mcmc, traceLog, new ArrayList<>(), null));
			table.getItems().add(new LoggerSelection(id, "Substitution count logger", substTreeLogger, mcmc, traceLog, loggers, substLogger));

			
		}
		
		
		
		
		// StarBeast3 ?
		SpeciesTree speciesTree = getSpeciesTree(likelihoods);
		if (speciesTree != null) {
			
			String id = speciesTree.getID();
			String substTreeLoggerID =  "BeastMap.SubstTreeLogger." + id;
			
			// Subst logger
			Logger substTreeLogger;
			//BranchMutationSampler sampler;
			if (inputContainsID(mcmc.loggersInput.get(), substTreeLoggerID)) {
			
				substTreeLogger = (Logger) getInputWithID(mcmc.loggersInput.get(), substTreeLoggerID);
				Log.warning("Found species logger in MCMC!");
			} else {
				Log.warning("Making new species logger called " + substTreeLoggerID);
			
				
				String outfile = "beastmap.species.trees";
				substTreeLogger = getTreeLoggerForLikelihood(speciesTree);
				
				
				// This will always be true unless there are morphological traits
				if (substTreeLogger == null) {
				
					// Loggable
					List<Loggable> substLoggables = new ArrayList<>();
					SampledSubstTreeLogger typedTreeLogger = new SampledSubstTreeLogger();
					typedTreeLogger.initByName("tree", speciesTree);
					typedTreeLogger.setID("BeastMap.SampledSubstTreeLogger." + id);
					substLoggables.add(typedTreeLogger);
					
					// Loggers
					substTreeLogger = new Logger();
					
					substTreeLogger.initByName("fileName", outfile, "logEvery", 10000, "mode", "tree", "log", substLoggables);
					
					
					doc.registerPlugin(typedTreeLogger);
					doc.registerPlugin(substTreeLogger);

					
				}
				
				
				
				// This logger should always be named after the species, not one of the morphological traits
				substTreeLogger.fileNameInput.set(outfile);
				substTreeLogger.setID(substTreeLoggerID);
				
				
			}
			
			SampledSubstTreeLogger substLogger = null;
			for (BEASTObject obj : substTreeLogger.loggersInput.get()) {
				if (obj instanceof SampledSubstTreeLogger) {
					substLogger = (SampledSubstTreeLogger) obj;
				}
			}
			
			
			
			
			List<SubstLoggerSelection> loggers = getSubstSumLoggersForSpeciesTree(speciesTree, id, substLogger, likelihoods, mappersAll);
						
			
			// Add this to the top of the table
			table.getItems().add(0, new LoggerSelection(id, "Substitution count logger", substTreeLogger, mcmc, traceLog, loggers, substLogger));
			
			
		}
		
		
		getChildren().add(pane);
    	
    }
    
    
    
    /**
     * See if there is an an existing beastmap logger with the same tree, and if not return null
     * @param likelihood
     * @return
     */
    private Logger getTreeLoggerForLikelihood(Tree tree ) {
    	
    	
    	
		for (String id : doc.pluginmap.keySet()) {
			
			BEASTInterface o = doc.pluginmap.get(id);
			
			// Is this a Logger?
			if (o instanceof Logger) {
				
				
				Logger logger = (Logger) o;
				
				// Does it have a SampledSubstTreeLogger?
				for (BEASTObject obj : logger.loggersInput.get()) {
					if (obj instanceof SampledSubstTreeLogger) {
						SampledSubstTreeLogger substLogger = (SampledSubstTreeLogger) obj;
						
						if (substLogger.likelihoodInput.get() == null) {
							Log.warning("Cannot find likelihood for " + substLogger.getID());
							continue;
						}
						
						// Same tree?
						if (substLogger.likelihoodInput.get().treeInput.get() == tree) {
							Log.warning("Found a preexisting tree logger for " + tree.getID());
							return logger;
						}
					}
				}
				
			}
			
		}
		
		return null;
    	
    	
    }
    
    

	private List<BranchMutationSampler> getAllMappers() {
		
		List<BranchMutationSampler> mappers = new ArrayList<>();
		Map<String, BEASTInterface> objects = doc.pluginmap;
		for (String key : objects.keySet()) {
			
			//Log.warning("Found " + key);
			BEASTInterface o = objects.get(key);
			
			// Is this a mapper?
			if (o instanceof BranchMutationSampler) {
				BranchMutationSampler b = (BranchMutationSampler)o;
				mappers.add(b);
				Log.warning("Found stochatstic mapper in mcmc " + b.getID());
			}
			
			
		}
		
		
		return mappers;
	}
	

	

	// Check if there is a species tree that embeds all gene trees (ie StarBeast3)
    private SpeciesTree getSpeciesTree(List<GenericTreeLikelihood> likelihoods) {
    	
    	//<distribution id="treePrior.t:10206_NT_AL" spec="starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution" populationModel="@speciesTreePopulationModel" speciesTree="@Tree.t:Species" speciesTreePrior="@SpeciesTreePopSize.Species" tree="@Tree.t:10206_NT_AL"/>
    	SpeciesTree speciesTree = null;
    	for (GenericTreeLikelihood likelihood : likelihoods) {
    		
    		// Gene tree
    		Tree tree = (Tree) likelihood.treeInput.get();
    		
    		// Gene tree prior
    		for (BEASTInterface obj : tree.getOutputs()) {
    			
    			
    			// Get a species tree (and check all genes have the same species tree)
    			if (obj instanceof GeneTreeForSpeciesTreeDistribution) {
    				
    				GeneTreeForSpeciesTreeDistribution prior = (GeneTreeForSpeciesTreeDistribution) obj;
    				SpeciesTree speciesTree2 = prior.speciesTreeInput.get();
    				
    				if (speciesTree == null) {
    					speciesTree = speciesTree2;
    				}
    				
    				if (speciesTree != speciesTree2) return null;
    				
    				
    			}
    			
    		}
    		
    		
    	}
    	
    	
    	if (speciesTree != null) {
    		Log.warning("Found species tree " + speciesTree.getID());
    	}
    	
    	return speciesTree;
    	
    }

    
    

    private List<SubstLoggerSelection> getSubstSumLoggers(BranchMutationSampler sampler, String id, SampledSubstTreeLogger logger, DataType dt)  {
    	
    	List<SubstLoggerSelection> loggers = new ArrayList<>();
    	for (Class<?> clazz : termsToLog) {
    		
    		
    		try {
				Object instance = clazz.getDeclaredConstructor().newInstance();
				BranchSubstLogger sum = (BranchSubstLogger) instance;
				
				if (!sum.canHandleDataType(dt)) {
					continue;
				}
				
				sum.initByName("sampler", sampler);
				sum.setID("beastmap." + sum.getClass().getSimpleName() + "." + id);
				SubstLoggerSelection obj = new SubstLoggerSelection(sum.getClass().getSimpleName(), true, sum);
				
				

	    		// See if the object already exists
	    		for (StochasticMapProperty term : logger.samplerInput.get()) {
	    			
	    			if (term instanceof BranchSubstLogger) {
	    				
	    				BranchSubstLogger term2 = (BranchSubstLogger) term;
	    				if (term2.getID().equals(sum.getID())) {
	    					obj = new SubstLoggerSelection(term.getClass().getSimpleName(), true, term2);
	    					Log.warning("Found the mut " + term2.getID());
	    					break;
	    				}
	    				
	    			}
	    			
	    		}
	    		
    			loggers.add(obj);
	    		
				
			} catch (Exception e) {
				
				e.printStackTrace();
				
			}
    		
    	}
    	
		return loggers;
    	
	}
    
    
    
    
    // Get all loggers for the species tree (summing over all gene trees)
    private List<SubstLoggerSelection> getSubstSumLoggersForSpeciesTree(SpeciesTree speciesTree, String id, SampledSubstTreeLogger substLogger, List<GenericTreeLikelihood> likelihoods, List<BranchMutationSampler> mappersAll) {
    	
    	

    	List<SubstLoggerSelection> loggers = new ArrayList<>();
    	for (Class<?> clazz : termsToLog) {
    		
    		
    		try {
    			
    			SubstitutionSummerPerBranch speciesSummer = new SubstitutionSummerPerBranch();
				String className = clazz.getDeclaredConstructor().newInstance().getClass().getSimpleName();
				speciesSummer.setID("beastmap.species." + className);
				
				
				List<BranchSubstLogger> summerList = new ArrayList<>();
				
				// Which gene trees can count this term?
				for (GenericTreeLikelihood likelihood : likelihoods) {
					
					
					
					Object instance = clazz.getDeclaredConstructor().newInstance();
					BranchSubstLogger sumGeneTree = (BranchSubstLogger) instance;
					
		    		// Gene tree
		    		Tree tree = (Tree) likelihood.treeInput.get();
		    		
		    		
		    		// Morphology
		    		if (tree == speciesTree) continue;
		    		
		    		// Valid data type?
		    		DataType dt = likelihood.dataInput.get().getDataType();
		    		if (!sumGeneTree.canHandleDataType(dt)) {
		    			//Log.warning(instance.getClass().getSimpleName() + " has the wrong data type " + dt.getClass().getSimpleName());
						continue;
					}
					
		    		
		    		// Sampler
		    		String samplerID = "BeastMap.StochasticMapper." + BeautiDoc.parsePartition(likelihood.getID());
		    		BranchMutationSampler sampler = (BranchMutationSampler) getInputWithID(mappersAll, samplerID);
		    		if (sampler == null) {
		    			Log.warning("Cannot find sampler with ID " + samplerID );
		    		}
		    		
		    		
		    		// Gene tree prior
		    		GeneTreeForSpeciesTreeDistribution prior = null;
		    		for (BEASTInterface obj : tree.getOutputs()) {
		    			if (obj instanceof GeneTreeForSpeciesTreeDistribution) {
		    				prior = (GeneTreeForSpeciesTreeDistribution) obj;
		    				SpeciesTree speciesTree2 = prior.speciesTreeInput.get();
		    				if (speciesTree != speciesTree2) {
		    					Log.warning("Unexpected error: " + speciesTree.getID() + " != " + speciesTree2);
		    				}
		    				break;
		    			}
		    		}
		    		if (prior == null) continue; // Don't know how this is possible...
		    		
		    		
		    		// Add the counter to the species summer list
		    		sumGeneTree.initByName("sampler", sampler, "includeRoot", true, "gene", prior);
		    		sumGeneTree.setID("beastmap.gene." + sumGeneTree.getClass().getSimpleName() + "." + likelihood.getID());
		    		
		    		
		    		//Log.warning(" sumGeneTree " + likelihood.getID() + " " + sumGeneTree.getID());
		    		
		    		summerList.add(sumGeneTree);
		    		
				}
				
				if (summerList.isEmpty()) {
					Log.warning(className + " does not have any summers");
					continue;
				}
				
				
				speciesSummer.initByName("counter", summerList);
				
				
				SubstLoggerSelection obj = new SubstLoggerSelection("Species " + className, true, speciesSummer);
				
				
				// See if the object already exists
	    		for (StochasticMapProperty term : substLogger.samplerInput.get()) {
	    			
	    			if (term instanceof SubstitutionSummerPerBranch) {
	    				
	    				SubstitutionSummerPerBranch term2 = (SubstitutionSummerPerBranch) term;
	    				if (term2.getID().equals(speciesSummer.getID())) {
	    					obj = new SubstLoggerSelection("Species " + className, true, term2);
	    					speciesSummer = term2;
	    					Log.warning("Found the species mut " + term2.getID());
	    					break;
	    				}
	    				
	    			}
	    			
	    		}
	    		
	    		
				
				loggers.add(obj);
				
	    		
				
			} catch (Exception e) {
				
				e.printStackTrace();
				
			}
    		
    	}
    	
		return loggers;
    	
		
	}
    

    
    
    private HBox prepareBurninTextbox(List<BranchMutationSampler> mappersAll) {
    	
    	
    	HBox row = new HBox(10);   // 10px spacing between items
    	
    	
    	Label nameLabel = new Label("Map burnin:");
    	TextField tf = new TextField();
    	
    	nameLabel.setPrefWidth(100); 
    	tf.setPrefWidth(100);
    	Tooltip tip = new Tooltip("The number of MCMC states until stochastic mapping is used. The algorithm is prone to numerical errors during the start of a chain, as the branch lengths can be very long. During burnin, there will be 0 substitutions assinged to every branch.");
    	tf.setTooltip(tip);
    	
    	row.getChildren().addAll(nameLabel, tf);
    	//row.setPadding(new Insets(0, 0, 0, 50)); // top, right, bottom, left
    	
    	
    	// Load initial value
    	int initialValue = burnin.get();
    	for (BranchMutationSampler sampler : mappersAll) {
    		initialValue = sampler.burninInput.get();
    		Log.warning("Loading burnin from " + sampler.getID() + " " + initialValue);
    		break;
    	}
    	
    	tf.setText("" + initialValue);
    	

		tf.textProperty().addListener((obs, old, nw) -> {
		    if (!nw.matches("\\d*(\\.\\d+)?")) {
		        tf.setText(old);
		    } else if (!nw.isEmpty()) {
		    	
		    	int val = Integer.parseInt(nw);
		    	val = Math.max(0, val);
		    	tf.setText("" + val);
		        burnin.set(val);
		        
		        for (BranchMutationSampler sampler : mappersAll) {
		        	
		        	sampler.burninInput.set(val);
		        	Log.warning("Updating burnin in " + sampler.getID() + " to " + sampler.burninInput.get());
		        }
		        
		        
		    }
		});
		
		return row;

    	
    }
    
    
	// Prepare the outer table
    private TableView<LoggerSelection> prepareTable(){
    	

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
	    	
    	
	    
	    // Name column
	    TableColumn<LoggerSelection, String> nameCol = new TableColumn<>("Partition");
	    nameCol.setPrefWidth(250);
	    nameCol.setCellValueFactory(data -> data.getValue().nameProperty());
	    table.getColumns().add(nameCol);
	    
	    
	    // Logger column
	    TableColumn<LoggerSelection, String> loggerCol = new TableColumn<>("Logger type");
	    loggerCol.setPrefWidth(250);
	    loggerCol.setCellValueFactory(data -> data.getValue().loggerProperty());
	    table.getColumns().add(loggerCol);
	    
	    
	    // Check boxes
	    TableColumn<LoggerSelection, Boolean> selectedCol = new TableColumn<>("Tree logger");
	    selectedCol.setGraphic(selecteAllCheckBox);
	    selectedCol.setSortable(false);
	    selectedCol.setPrefWidth(130);
	    

	    selectedCol.setCellFactory( tc -> {
	    	CheckBoxTableCell cell =  new CheckBoxTableCell<>();
	    	cell.setEditable(true);
	    	List o = cell.getChildrenUnmodifiable();
	    	return cell;
	    });
	    selectedCol.setCellValueFactory( f -> f.getValue().selectedProperty());
	    
	    selectedCol.setCellFactory(CheckBoxTableCell.forTableColumn(selectedCol));

	    
	    selectedCol.setEditable(true);
	    table.getColumns().add(selectedCol);
	    
	   

	    
	    TableColumn<LoggerSelection, List<String>> optionsColumn = new TableColumn<>("Terms to count");
	    optionsColumn.setPrefWidth(500);
	    optionsColumn.setCellValueFactory(data ->
		    new ReadOnlyObjectWrapper<>(data.getValue().getOptionNames())
		);
	    
	    
	    
	    optionsColumn.setCellFactory(col -> new TableCell<LoggerSelection, List<String>>() {

	        private final VBox box = new VBox(4); // spacing between checkboxes

	        @Override
	        protected void updateItem(List<String> items, boolean empty) {
	            super.updateItem(items, empty);

	            if (empty || items == null) {
	                setGraphic(null);
	            } else {
	                box.getChildren().clear();

	                for (String name : items) {
	                    CheckBox cb = new CheckBox(name);
	                    
	                    
	                    // Bind check state to the row’s underlying data:
	                    LoggerSelection row = getTableView().getItems().get(getIndex());

	                    
	                    cb.setSelected(row.isChecked(name));
	                    cb.selectedProperty().addListener((obs, oldV, newV) ->
	                        row.setChecked(name, newV)
	                    );
	                    
	                    

	                    box.getChildren().add(cb);
	                    
	                }

	                setGraphic(box);
	            }
	        }
	    });
	    
	    table.getColumns().add(optionsColumn);
	    
	    
	    return table;
    	
    }
    

    
    public void refreshPanel() {
       
    }
    
    public static boolean customConnector(BeautiDoc doc) {
    	
    	
    	System.out.println("BEASTMAP HERE");
    	return true;
    	
    }
    
    

	// Add to list if not already there
	public static <T> BEASTObject addWithID(List<T> list, BEASTObject toAdd) {
		
		for (T obj : list) {
			
			if (obj instanceof BEASTObject) {
				
				BEASTObject bobj = (BEASTObject) obj;
				if (bobj.getID().equals(toAdd.getID())) {
					return bobj;
				}
			}
			
		}
		
		list.add((T) toAdd);
		
		
		return toAdd;
		
		
	}
	
	
	
	// Remove the input if it has the right id
	public <T> void removeWithID(List<T> list, String id) {
		
		
		
		Iterator<T> iter = list.iterator();
		while (iter.hasNext()) {
			
			T obj = iter.next();
			
			if (obj instanceof BEASTObject) {
				
				BEASTObject bobj = (BEASTObject) obj;
				if (bobj.getID().equals(id)) {
					iter.remove();
					//doc.unregisterPlugin(bobj);
				}
			}
			
		}
		

		
	}
	
	
	// Does the input contain something with this ID?
	public <T> boolean inputContainsID(List<T> list, String id) {
		
		for (T obj : list) {
			
			if (obj instanceof BEASTObject) {
				
				BEASTObject bobj = (BEASTObject) obj;
				if (bobj.getID().equals(id)) {
					return true;
				}
			}
			
		}
		
		return false;
		
	}
	
	
	// Does the input contain something with this ID?
	public <T> BEASTObject getInputWithID(List<T> list, String id) {
		
		for (T obj : list) {
			
			if (obj instanceof BEASTObject) {
				
				BEASTObject bobj = (BEASTObject) obj;
				if (bobj.getID().equals(id)) {
					return bobj;
				}
			}
			
		}
		
		return null;
		
	}
    
    

    // Include a logger in MCMC?
    public class LoggerSelection {
    	
		private StringProperty name = new SimpleStringProperty();
		private StringProperty loggerType = new SimpleStringProperty();
		private MCMC mcmc;
		private Logger traceLog;
		private Logger logger;
	    private BooleanProperty selected = new SimpleBooleanProperty();
	    private List<SubstLoggerSelection> options;
	    private SampledSubstTreeLogger substLogger;
	    
	    
	    public LoggerSelection(String name, String loggerType, Logger logger, MCMC mcmc, Logger traceLog, List<SubstLoggerSelection> options, SampledSubstTreeLogger substLogger) {
	    	
	    	this.setName(name);
	    	this.setLoggerType(loggerType);
	    	this.logger = logger;
	    	this.mcmc = mcmc;
	    	this.traceLog = traceLog;
	    	this.options = options;
	    	this.substLogger = substLogger;
	    	
	    	
	    	setSelected(inputContainsID(this.mcmc.loggersInput.get(), logger.getID()));
	    	selected.addListener((obs,  wasSelected,  isSelected) -> {
				setSelected(isSelected);
			});	
			
			
			// Which terms are being logged?
			for (SubstLoggerSelection option : this.options) {
				boolean active = inputContainsID(substLogger.samplerInput.get(), option.getLoggerID()) || doc.pluginmap.containsKey(option.getLoggerID());
				Log.warning("Setting " + option.getName() + "to " + active);
				this.setChecked(option.getName(), active);
			}
			
			
	    }
	    
	    


		public Logger getLogger() {
			return logger;
		}



		// Is the logger on?
	    public Boolean isSelected() { return selected.get(); }
	    public BooleanProperty selectedProperty() { return selected; }
	    
	    
	    
	    // Same as setSelected but without any recursion
	    public void resetSelection() {
	    	
			boolean s = inputContainsID(this.mcmc.loggersInput.get(), logger.getID());
			this.selected.set(s); 
			
			if (s) {
	    		addWithID(this.mcmc.loggersInput.get(), logger);
	    	}else {
	    		removeWithID(this.mcmc.loggersInput.get(), logger.getID());
	    	}
			
	    }
	    
	    public void setSelected(Boolean selected) { 
	    	
	    	
	    	//Log.warning("SELECTED " + this.logger.getID());
	    	
	    	this.selected.set(selected); 
			
	    	if (selected) {
	    		addWithID(this.mcmc.loggersInput.get(), logger);
	    	}else {
	    		removeWithID(this.mcmc.loggersInput.get(), logger.getID());
	    	}
	    	
	    	
	    	for (LoggerSelection item : table.getItems()) {
	    		item.resetSelection();
	    	}
	    	
	    }
	    

	    
	    // What terms are being logged by the subst logger?
	    public List<String> getOptionNames() {
			List<String> opts = new ArrayList<>();
			for (SubstLoggerSelection option : this.options) {
				opts.add(option.getName());
			}
			return opts;
		}

		public void setChecked(String name, Boolean value) {
			
			for (SubstLoggerSelection option : this.options) {
				if (option.getName().equals(name)) {
					option.setActive(value);
					
					
					// Add/remove the term to the tree logger
					if (value) {
						
						addWithID(substLogger.samplerInput.get(), (BEASTObject)option.loggable);
						
//						if (!substLogger.samplerInput.get().contains(option.logger)) {
//							substLogger.samplerInput.get().add(option.logger);
//						}
						
						doc.registerPlugin((BEASTObject)option.loggable);
						
					}else {
						//substLogger.samplerInput.get().remove(option.logger);
						
						BEASTObject o = getInputWithID(substLogger.samplerInput.get(), option.getLoggerID());
						removeWithID(substLogger.samplerInput.get(), option.getLoggerID());
						
						if (o != null){
							doc.unregisterPlugin((BEASTObject)o);
						}
						
						
					}
					
					
					// Add/remove their sums to the trace logger (if appropriate)
					SubstitutionSummer sum_summer = new SubstitutionSummer();
					sum_summer.initByName("counter", option.loggable);
					String newID = "beastmap.sumof." + option.name + "." + this.getName();
					sum_summer.setID(newID);
					if (value) {
						
						boolean alreadyThere = false;
						for (BEASTObject l : traceLog.loggersInput.get()) {
							if (l.getID().equals(newID)) {
								alreadyThere = true;
								break;
							}
						}
						
						if (!alreadyThere) {
							traceLog.loggersInput.get().add(sum_summer);
						}
						
						
						
					}else {
						
						for (BEASTObject l : traceLog.loggersInput.get()) {
							if (l.getID().equals(newID)) {
								traceLog.loggersInput.get().remove(l);
								break;
							}
						}
						
					}
					
					return;
				}
			}
		}
		
		
		

		
		// Are we counting this term? eg SubstitutionSum
		public boolean isChecked(String name) {
			
			//if (!isSelected()) return false;
			
			for (SubstLoggerSelection option : this.options) {
				if (option.getName().equals(name)) {
					return option.isActive();
				}
			}
			
			return false;
		}


		
	   
	    public void setName(String name) { this.name.set(name); }
	    public String getName() { return name.get(); }
	    public StringProperty nameProperty() { return name; }
	    
	    public void setLoggerType(String type) { this.loggerType.set(type); }
	    public StringProperty loggerProperty() { return loggerType; }
	
	   
	    
    }
    
    
    
    

    // Include a logger in MCMC?
    public class SubstLoggerSelection {
    	
    	String name;
    	boolean active;
    	Loggable loggable;
    	
    	public SubstLoggerSelection(String name, boolean active, Loggable logger) {
    		this.name = name.replaceAll(" ", "");
    		this.active = active;
    		this.loggable = logger;
    	}



		public String getName() {
			return name;
		}
		
		
		public void setActive(Boolean value) {
			this.active = value;
		}
		
		public boolean isActive() {
			return active;
		}
		
		public Loggable getLogger() {
			return loggable;
		}
		
		
		public String getLoggerID() {
			return ((BEASTObject) loggable).getID();
		}
    	
	    
    }
    
    
    
    
}


