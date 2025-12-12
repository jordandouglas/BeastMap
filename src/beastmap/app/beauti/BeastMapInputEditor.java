package beastmap.app.beauti;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
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
import beast.base.core.Loggable;
import beastmap.logger.SampledSubstTreeLogger;
import beastmap.logger.StochasticMapProperty;
import beastmap.logger.SubstitutionSummer;






public class BeastMapInputEditor extends ListInputEditor {

	
	public static Class[] termsToLog = new Class[] { AminoAcidClassChanges.class, AminoAcidClassRemains.class, 
										 NonSynonymousSubstSum.class, NucleotideTransitionCounter.class, 
										NucleotideTransversionCounter.class, SubstitutionSum.class, SynonymousSubstSum.class,  }; // TotalSize.class NetSize.class, FromToSubstSum.class,
	
	
	
	static HashMap<String, BranchMutationSampler> mapperMap = new HashMap<>();
	static HashMap<String, Logger> segmentedMap = new HashMap<>();
	static HashMap<String, Logger> substMap = new HashMap<>();
	
	
	private final SimpleIntegerProperty burnin = new SimpleIntegerProperty(100000);
	
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
		label2.getChildren().add(new Label("\t•Segmented tree loggers produce one branch segment each time the state changes (ideal for discrete geography). Warning: tree files may grow very large."));
		HBox label3 = FXUtils.newHBox();
		label3.getChildren().add(new Label("\t•Branch substitution loggers count the number of events on each branch."));
		pane.getChildren().add(label1);
		pane.getChildren().add(label2);
		pane.getChildren().add(label3);
		
		
		pane.getChildren().add(prepareBurninTextbox());
		
		
	    // Add table to page
		TableView<LoggerSelection> table = this.prepareTable();
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
				
		
		
		// Remove any unused likelihoods from the cache
		for (String id : mapperMap.keySet()) {
			
			Log.warning("id=" + id);
			
			boolean likelihoodActive = false;
			for (GenericTreeLikelihood likelihood : likelihoods) {
				if (likelihood.getID().equals(id)) {
					likelihoodActive = true;
					break;
				}
			}

			if (!likelihoodActive) {
				Log.warning("Removing " + id);
				
				for (Logger logger : mcmc.loggersInput.get()) {
					if (logger.getID().equals(segmentedMap.get(id).getID())) {
						mcmc.loggersInput.get().remove(logger);
					}else if (logger.getID().equals(substMap.get(id).getID())) {
						mcmc.loggersInput.get().remove(logger);
					}
				}
				
				segmentedMap.remove(id);
				substMap.remove(id);
				mapperMap.remove(id);
				
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
			
				Alignment data = likelihood.dataInput.get();
				
				
				// Dataset
				PatternlessAlignment dataPatternless = new PatternlessAlignment();
				dataPatternless.initByName("data", data); // User data type??;
				dataPatternless.setID("BeastMap.PatternlessAlignment." + id);
				
				
				// Sampler
				sampler = new BranchMutationSampler();
				sampler.initByName("likelihood", likelihood, "burnin", 100000, "data", dataPatternless);
				sampler.setID("BeastMap.BranchMutationSampler." + id);
				
				
				
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
				segmentedTreeLogger.setID("BeastMap.SegmentedTreeLogger." + id);
				
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

				
				// Samplers
//				List<BranchSubstLogger> substCounters = new ArrayList<>();
//				SubstitutionSum sum = new SubstitutionSum();
//				sum.initByName("sampler", sampler);
//				sum.setID("BeastMap.SubstitutionSum." + id);
//				substCounters.add(sum);
				
				
				
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
				substTreeLogger.setID("BeastMap.SubstTreeLogger." + id);
				
				// Add to MCMC chain
				//mcmc.loggersInput.get().add(segmentedTreeLogger);
				substMap.put(likelihood.getID(), substTreeLogger);
				
				

				
				Log.warning("Adding logger to map!");
				
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
		
		
		getChildren().add(pane);
    	
    }
    

    
    

    private List<SubstLoggerSelection> getSubstSumLoggers(BranchMutationSampler sampler, String id, SampledSubstTreeLogger logger, DataType dt)  {
    	
    	
    	
    	
    	

    	List<SubstLoggerSelection> loggers = new ArrayList<>();
    	for (Class<?> clazz : termsToLog) {
    		
    		
    		try {
				Object instance = clazz.getDeclaredConstructor().newInstance();
				BranchSubstLogger sum = (BranchSubstLogger) instance;
				sum.initByName("sampler", sampler);
				sum.setID("beastmap." + sum.getClass().getSimpleName() + "." + id);
				SubstLoggerSelection obj = new SubstLoggerSelection(sum.getClass().getSimpleName(), true, sum);
				
				

	    		// See if the object already exists
	    		for (StochasticMapProperty term : logger.samplerInput.get()) {
	    			
	    			if (term instanceof BranchSubstLogger) {
	    				
	    				BranchSubstLogger term2 = (BranchSubstLogger) term;
	    				if (term2.getID().equals(sum.getID())) {
	    					obj = new SubstLoggerSelection(term.getClass().getSimpleName(), true, term2);
	    					break;
	    				}
	    				
	    			}
	    			
	    		}
	    		
	    		// Does it handle this datatype?
	    		if (sum.canHandleDataType(dt)) {
	    			loggers.add(obj);
	    		}
	    		
				
			} catch (Exception e) {
				
				e.printStackTrace();
				
			}
    		

    		
    	}
    	
    
		
		
		return loggers;
    	
	}

    
    
    private HBox prepareBurninTextbox() {
    	
    	
    	HBox row = new HBox(10);   // 10px spacing between items
    	
    	
    	Label nameLabel = new Label("Map burnin: ");
    	TextField tf = new TextField();
    	
    	nameLabel.setPrefWidth(100); 
    	tf.setPrefWidth(100);
    	Tooltip tip = new Tooltip("The number of MCMC states until stochastic mapping is used. The algorithm is prone to numerical errors during the start of a chain, as the branch lengths can be very long. During burnin, there will be 0 substitutions assinged to every branch.");
    	tf.setTooltip(tip);
    	
    	row.getChildren().addAll(nameLabel, tf);
    	//row.setPadding(new Insets(0, 0, 0, 50)); // top, right, bottom, left
    	
    	
    	// Load initial value
    	int initialValue = burnin.get();
    	for (String id : mapperMap.keySet()) {
    		initialValue = mapperMap.get(id).burninInput.get();
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
		        
		        for (String id : mapperMap.keySet()) {
		        	mapperMap.get(id).burninInput.set(val);
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
	    TableColumn<LoggerSelection, Boolean> selectedCol = new TableColumn<>("Active");
	    selectedCol.setGraphic(selecteAllCheckBox);
	    selectedCol.setSortable(false);
	    selectedCol.setPrefWidth(100);
	    

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


	    
	    // Logger settings column
//	    TableColumn<LoggerSelection, SubstLoggerSelection> settingsCol = new TableColumn<>("Terms to count");
//	    settingsCol.setPrefWidth(500);
//	    settingsCol.setCellValueFactory(data -> data.getValue().substLoggerTable());
//	    table.getColumns().add(settingsCol);
	    
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
	                    
	                    

	                    // optional: bind check state to the row’s underlying data:
	                    LoggerSelection row = getTableView().getItems().get(getIndex());
	    
                    	// replace row.isChecked(name) / row.setChecked(name, value) with your getters/setters
	                    cb.setSelected(row.isChecked(name));
	                    cb.selectedProperty().addListener((obs, oldV, newV) ->
	                        row.setChecked(name, newV)
	                    );
//	                    cb.disabledProperty().addListener((obs, oldV, newV) ->
//	                        row.isSelected()
//	                    );
//                    
	                    
	                 

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
	    	
	    	
	    	// Is the logger on?
			setSelected(this.mcmc.loggersInput.get().contains(logger));
			selected.addListener((obs,  wasSelected,  isSelected) -> {
				setSelected(isSelected);
			});	
			
			
			// Which terms are being logged?
			for (SubstLoggerSelection option : this.options) {
				boolean active = logger.loggersInput.get().contains(option.getLogger());
				this.setChecked(option.getName(), active);
			}
			
			
	    }
	    
	    
	    
	    // Is the logger on?
	    public Boolean isSelected() { return selected.get(); }
	    public BooleanProperty selectedProperty() { return selected; }
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
						if (!substLogger.samplerInput.get().contains(option.logger)) {
							substLogger.samplerInput.get().add(option.logger);
						}
						
					}else {
						substLogger.samplerInput.get().remove(option.logger);
					}
					
					
					// Add/remove their sums to the trace logger (if appropriate)
					SubstitutionSummer sum_summer = new SubstitutionSummer();
					sum_summer.initByName("counter", option.logger);
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
    	BranchSubstLogger logger;
    	
    	public SubstLoggerSelection(String name, boolean active, BranchSubstLogger logger) {
    		this.name = name;
    		this.active = active;
    		this.logger = logger;
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
		
		public BranchSubstLogger getLogger() {
			return logger;
		}
    	
		
	    
    }
    
    
    
    
}


