package beastmap.app.beauti;


import beast.base.core.Description;
import beast.base.core.Log;
import beast.base.core.ProgramStatus;

import java.io.File;
import java.util.*;

import beastfx.app.inputeditor.AlignmentImporter;
import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiSubTemplate;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import beastfx.app.util.Utils;
import beastmap.indel.SimpleIndelCodingAlignment;
import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.Binary;
import beast.base.parser.PartitionContext;
import beast.pkgmgmt.BEASTClassLoader;

@Description("Class for creating new partitions for simple indel representations, edited by AlignmentListInputEditor")
public class BeautiSimpleIndelAlignmentProvider extends BeautiAlignmentProvider {

	
	static List<AlignmentImporter> importers = null;
	
	
	private void initImporters() {
		importers = new ArrayList<>();		

        // build up list of data types
        Set<String> importerClasses = Utils.loadService(AlignmentImporter.class);        
        for (String _class: importerClasses) {
        	try {
        		if (!_class.startsWith(this.getClass().getName())) {
					AlignmentImporter importer = (AlignmentImporter) BEASTClassLoader.forName(_class).newInstance();
					importers.add(importer);
        		}
			} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }
	}
	
	

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		
		
		Log.warning("Get indel alignments!");
		
		if (files == null) {
			// merge "+ button" and "drag drop" function
			return getAlignments(doc);
		}
		if (importers == null) {
			initImporters();
		}
		
		
		List<BEASTInterface> selectedBEASTObjects = new ArrayList<>();
        for (File file : files) {
			// create list of importers that can handle the file
			List<AlignmentImporter> availableImporters = new ArrayList<>();
			for (AlignmentImporter importer : importers) {
				if (importer.canHandleFile(file)) {
					availableImporters.add(importer);
				}
			}
			if (availableImporters.size() == 0 && file.getPath().toLowerCase().endsWith(".txt")) {
				// remove .txt extension and try again
				String path = file.getPath();
				path = path.substring(0, path.length() - 4);
				File file2 = new File(path);
				for (AlignmentImporter importer : importers) {
					if (importer.canHandleFile(file2)) {
						availableImporters.add(importer);
					}
				}
			}
			
			if (availableImporters.size() > 0) {
				AlignmentImporter importer = availableImporters.get(0);
				if (availableImporters.size() > 1) {
					// let user choose an importer
					List<String> descriptions = new ArrayList<>();
					for (AlignmentImporter i : availableImporters) {
						descriptions.add(i.getDescription());
					}
					String option = (String)Alert.showInputDialog(null, "Which importer is appropriate", "Option",
		                    Alert.WARNING_MESSAGE, null, descriptions.toArray(), descriptions.get(0));
					if (option == null) {
						return selectedBEASTObjects;
					}
					int i = descriptions.indexOf(option);
					importer = availableImporters.get(i);
				}

				List<BEASTInterface> list = importer.loadFile(file);
				for (BEASTInterface o : list) {
					if (o.getID() != null && o.getID().contains(":")) {
						o.setID(o.getID().replaceAll(":", "-"));
					}
				}
				selectedBEASTObjects.addAll(list);
			} else {
                Alert.showMessageDialog(null,
                        "Unsupported sequence file.",
                        "Error", Alert.ERROR_MESSAGE);
			}
			
        }
		

        try {
        
			List<BEASTInterface> alignments = new ArrayList<>();
	        for (BEASTInterface o : selectedBEASTObjects) {
	            if (o instanceof Alignment) {
	            	
	            	// Add the original aligmment
	            	Alignment alignment = (Alignment) o;
	            	doc.addPlugin(alignment);
	            	
	            	BeautiSubTemplate originalTemplate = (BeautiSubTemplate) doc.pluginmap.get("StandardPartitionTemplate");
	            	String nameOriginal = alignment.getID();
	            	//PartitionContext contextOriginal = new PartitionContext(nameOriginal, nameOriginal, nameOriginal, nameOriginal);
	            	
	            	Log.warning(nameOriginal + "_" + originalTemplate.getID());
	            	doc.addAlignmentWithSubnet(alignment, originalTemplate);
	            	
	            	
	            	
	            	// Add the gapped alignment
	            	SimpleIndelCodingAlignment indelAlignment = new SimpleIndelCodingAlignment();
	            	indelAlignment.initByName("data", alignment, "statecount", 2, "datatype","binary");
	            	indelAlignment.setID(alignment.getID() + "_simpleindel");
	            	
	            	
	            	
	            	doc.addPlugin(indelAlignment);
	            	
	            	
	            	String nameIndel = indelAlignment.getID();
	            	PartitionContext contextIndel = new PartitionContext(nameIndel, nameIndel, nameIndel, nameOriginal);
	            	doc.addAlignmentWithSubnet(contextIndel, template.get());
	            	
	            	
	            	alignments.add(alignment);
	            	alignments.add(indelAlignment);
	            	
	            }
	        }
	        
	        Log.warning("returning alignments " + alignments.size());
		    return alignments;
        
        } catch (Exception e) {
        	e.printStackTrace();
        	return null;
        }
		
     
		
	}

	
	
	@Override
	public void editAlignment(Alignment alignment, BeautiDoc doc) {
		
		super.editAlignment(alignment, doc);
	}


	/** 
	 * return amount to which the provided matches an alignment 
	 * The provider with the highest match will be used to edit the alignment 
	 * */
	@Override
	public int matches(Alignment alignment) {
		
		if (alignment.getDataType() instanceof Binary) {
			return 20;
		}
		
		return 0;
	}



}
