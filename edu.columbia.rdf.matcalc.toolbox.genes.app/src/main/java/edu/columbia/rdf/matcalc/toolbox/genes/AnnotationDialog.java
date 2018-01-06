package edu.columbia.rdf.matcalc.toolbox.genes;

import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.Box;

import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.text.Splitter;
import org.jebtk.modern.UI;
import org.jebtk.modern.button.ModernCheckSwitch;
import org.jebtk.modern.dialog.ModernDialogHelpWindow;
import org.jebtk.modern.dialog.ModernDialogStatus;
import org.jebtk.modern.dialog.ModernMessageDialog;
import org.jebtk.modern.event.ModernClickEvent;
import org.jebtk.modern.panel.HBox;
import org.jebtk.modern.panel.HExpandBox;
import org.jebtk.modern.panel.ModernPanel;
import org.jebtk.modern.panel.VBox;
import org.jebtk.modern.spinner.ModernCompactSpinner;
import org.jebtk.modern.text.ModernAutoSizeLabel;
import org.jebtk.modern.text.ModernClipboardTextField;
import org.jebtk.modern.text.ModernTextBorderPanel;
import org.jebtk.modern.text.SuggestionTextBox;
import org.jebtk.modern.widget.ModernTwoStateWidget;
import org.jebtk.modern.widget.ModernWidget;
import org.jebtk.modern.window.ModernWindow;
import org.jebtk.modern.window.WindowWidgetFocusEvents;

import edu.columbia.rdf.matcalc.bio.AnnotationSidePanel;

/**
 * Control which conservation scores are shown.
 * 
 * @author Antony Holmes Holmes
 *
 */
public class AnnotationDialog extends ModernDialogHelpWindow {
  private static final long serialVersionUID = 1L;

  // private ModernRadioButton mCheckHuman =
  // new ModernRadioButton("Human");

  // private ModernRadioButton mCheckMouse =
  // new ModernRadioButton("Mouse");

  private ModernTwoStateWidget mCheckOverlappingGenes = new ModernCheckSwitch("Overlapping", true);

  private ModernTwoStateWidget mCheckClosestGene = new ModernCheckSwitch("1st closest", true);

  private ModernTwoStateWidget mCheck2ndClosestGene = new ModernCheckSwitch("2nd closest");

  private ModernTwoStateWidget mCheck3rdClosestGene = new ModernCheckSwitch("3rd closest");

  private ModernTwoStateWidget mCheck4thClosestGene = new ModernCheckSwitch("4th closest");

  private ModernTwoStateWidget mCheck5thClosestGene = new ModernCheckSwitch("5th closest");

  private ModernTwoStateWidget mCheckOtherClosest = new ModernCheckSwitch("nth closest");

  private ModernCompactSpinner mTextExt5p = new ModernCompactSpinner(0, 100000, 2000, 1000, false);

  private ModernCompactSpinner mTextExt3p = new ModernCompactSpinner(0, 100000, 1000, 1000, false);

  private ModernClipboardTextField mTextOtherClosest = new SuggestionTextBox("e.g. 1,2,4-6");

  // private SpeciesCombo mSpeciesCombo;

  private AnnotationSidePanel mGenomesPanel = new AnnotationSidePanel();

  public AnnotationDialog(ModernWindow parent) {
    super(parent, "geneannotation.help.url");

    setTitle("Gene Annotation");

    createUi();

    setup();
  }

  private void setup() {
    setResizable(true);

    mTextOtherClosest.addFocusListener(new FocusListener() {

      @Override
      public void focusGained(FocusEvent arg0) {
        mCheckOtherClosest.setSelected(true);
      }

      @Override
      public void focusLost(FocusEvent arg0) {
        // TODO Auto-generated method stub

      }
    });

    addWindowListener(new WindowWidgetFocusEvents(mOkButton));

    setSize(700, 480);

    UI.centerWindowToScreen(this);
  }

  private final void createUi() {
    // this.getContentPane().add(new JLabel("Change " +
    // getProductDetails().getProductName() + " settings", JLabel.LEFT),
    // BorderLayout.PAGE_START);

    Box box2;

    Box box = VBox.create();

    box2 = HBox.create();
    box2.add(new ModernAutoSizeLabel("5p"));
    box2.add(ModernPanel.createHGap());
    box2.add(mTextExt5p);
    // box2.add(ModernPanel.createHGap());
    // box2.add(new ModernAutoSizeLabel("bp"));
    box2.add(UI.createHGap(20));
    box2.add(new ModernAutoSizeLabel("3p"));
    box2.add(ModernPanel.createHGap());
    box2.add(mTextExt3p);
    box2.add(ModernPanel.createHGap());
    box2.add(new ModernAutoSizeLabel("bp"));
    box.add(new HExpandBox("Promoter", box2));

    box.add(UI.createVGap(10));

    // midSectionHeader("Genes", box);
    box.add(mCheckOverlappingGenes);
    // box.add(UI.createVGap(10));
    // box.add(UI.createVGap(5));
    box.add(mCheckClosestGene);
    // box.add(UI.createVGap(5));
    box.add(mCheck2ndClosestGene);
    // box.add(UI.createVGap(5));
    box.add(mCheck3rdClosestGene);
    // box.add(UI.createVGap(5));
    box.add(new HExpandBox(mCheckOtherClosest, new ModernTextBorderPanel(mTextOtherClosest, 250)));
    // box.add(UI.createVGap(5));
    // box.add(mCheck4thClosestGene);
    // box.add(UI.createVGap(5));
    // box.add(mCheck5thClosestGene);
    // box.add(UI.createVGap(5));
    box.setBorder(ModernWidget.BORDER);

    setCard(box);

    getTabsPane().addLeftTab("Genomes", mGenomesPanel, 200, 100, 400);
  }

  @Override
  public void clicked(ModernClickEvent e) {
    if (e.getSource().equals(mOkButton)) {
      if (getGenomes().size() == 0) {
        ModernMessageDialog.createWarningDialog(mParent, "You must select at least one annotation set.");
        return;
      }

    }

    super.clicked(e);
  }

  public boolean getShowOverlappingGenes() {
    return mCheckOverlappingGenes.isSelected();
  }

  public boolean getShowClosestGene() {
    return mCheckClosestGene.isSelected();
  }

  public boolean getShow2ndClosestGene() {
    return mCheck2ndClosestGene.isSelected();
  }

  public boolean getShow3rdClosestGene() {
    return mCheck3rdClosestGene.isSelected();
  }

  public boolean getShow4thClosestGene() {
    return mCheck4thClosestGene.isSelected();
  }

  public boolean getShow5thClosestGene() {
    return mCheck5thClosestGene.isSelected();
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

  public List<Integer> getClosestList() {
    Set<Integer> ret = new HashSet<Integer>();

    if (mCheckClosestGene.isSelected()) {
      ret.add(1);
    }

    if (mCheck2ndClosestGene.isSelected()) {
      ret.add(2);
    }

    if (mCheck3rdClosestGene.isSelected()) {
      ret.add(3);
    }

    if (mCheckOtherClosest.isSelected()) {
      for (String term : Splitter.onComma().text(mTextOtherClosest.getText())) {
        if (term.contains("-")) {
          List<String> tokens = Splitter.on("-").text(term);

          int s = Integer.parseInt(tokens.get(0));
          int e = Integer.parseInt(tokens.get(1));

          if (s > e) {
            int t = e;

            e = s;
            s = t;
          }

          // Can't use numbers below 1
          s = Math.max(1, s);

          for (int i = s; i <= e; ++i) {
            ret.add(i);
          }
        } else {
          ret.add(Math.max(1, Integer.parseInt(term)));
        }
      }
    }

    return CollectionUtils.sort(ret);
  }
}
