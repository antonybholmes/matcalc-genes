package edu.columbia.rdf.matcalc.toolbox.genes;

import java.util.List;

import javax.swing.Box;

import org.jebtk.modern.UI;
import org.jebtk.modern.dialog.ModernDialogHelpWindow;
import org.jebtk.modern.panel.HBox;
import org.jebtk.modern.panel.HExpandBox;
import org.jebtk.modern.panel.ModernPanel;
import org.jebtk.modern.panel.VBox;
import org.jebtk.modern.spinner.ModernCompactSpinner;
import org.jebtk.modern.text.ModernAutoSizeLabel;
import org.jebtk.modern.window.ModernWindow;
import org.jebtk.modern.window.WindowWidgetFocusEvents;

import edu.columbia.rdf.matcalc.bio.GenomesPanel;

/**
 * Control which conservation scores are shown.
 * 
 * @author Antony Holmes Holmes
 *
 */
public class TSSAnnotationDialog extends ModernDialogHelpWindow {
	private static final long serialVersionUID = 1L;
	
	private ModernCompactSpinner mTextExt5p = 
			new ModernCompactSpinner(0, 100000, 2000, 1000, false);
	
	private ModernCompactSpinner mTextExt3p = 
			new ModernCompactSpinner(0, 100000, 2000, 1000, false);

	//private SpeciesCombo mSpeciesCombo;


	private GenomesPanel mGenomesPanel = new GenomesPanel();
	
	public TSSAnnotationDialog(ModernWindow parent) {
		super(parent, "geneannotation.help.url");
		
		setTitle("TSS Annotation");

		createUi();
		
		setup();
	}

	private void setup() {
		addWindowListener(new WindowWidgetFocusEvents(mOkButton));
		
		setSize(640, 500);
		
		UI.centerWindowToScreen(this);
	}

	private final void createUi() {
		//this.getContentPane().add(new JLabel("Change " + getProductDetails().getProductName() + " settings", JLabel.LEFT), BorderLayout.PAGE_START);

		Box box2;
		
		Box box = VBox.create();
		sectionHeader("Genome", box);

		UI.setSize(mGenomesPanel, 600, 240);
		box.add(mGenomesPanel);
		
		midSectionHeader("Buffer", box);
		
		box2 = HBox.create();
		box2.add(new ModernAutoSizeLabel("5p"));
		box2.add(ModernPanel.createHGap());
		box2.add(mTextExt5p);
		//box2.add(ModernPanel.createHGap());
		//box2.add(new ModernAutoSizeLabel("bp"));
		box2.add(UI.createHGap(20));
		box2.add(new ModernAutoSizeLabel("3p"));
		box2.add(ModernPanel.createHGap());
		box2.add(mTextExt3p);
		box2.add(ModernPanel.createHGap());
		box2.add(new ModernAutoSizeLabel("bp"));
		box.add(new HExpandBox("Promoter", box2));
		
		setCardContent(box);
	}

	public int getExt5p() {
		return Math.abs(mTextExt5p.getIntValue());
	}
	
	public int getExt3p() {
		return Math.abs(mTextExt3p.getIntValue());
	}

	public List<String> getGenomes() {
		return mGenomesPanel.getGenomes();
	}
}
