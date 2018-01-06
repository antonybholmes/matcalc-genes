package edu.columbia.rdf.matcalc.toolbox.genes;

import java.io.IOException;
import java.util.List;

import org.jebtk.modern.UIService;
import org.jebtk.modern.dialog.ModernDialogStatus;
import org.jebtk.modern.dialog.ModernMessageDialog;
import org.jebtk.modern.event.ModernClickEvent;
import org.jebtk.modern.event.ModernClickListener;
import org.jebtk.modern.graphics.icons.ModernIcon;
import org.jebtk.modern.help.ModernMenuHelpItem;
import org.jebtk.modern.menu.ModernPopupMenu;
import org.jebtk.modern.menu.ModernTwoLineMenuItem;
import org.jebtk.modern.ribbon.RibbonLargeDropDownButton;
import edu.columbia.rdf.matcalc.MainMatCalcWindow;
import edu.columbia.rdf.matcalc.toolbox.CalcModule;

public class GeneAnnotationModule extends CalcModule implements ModernClickListener {

  private static final ModernIcon ICON = UIService.getInstance().loadIcon("genes", 24);
  private MainMatCalcWindow mWindow;
  private RibbonLargeDropDownButton mButton;

  @Override
  public String getName() {
    return "Gene Annotation";
  }

  @Override
  public void init(MainMatCalcWindow window) {
    mWindow = window;

    ModernPopupMenu popup = new ModernPopupMenu();

    popup.addMenuItem(new ModernTwoLineMenuItem("Genes", "Annotate regions for overlapping genes.", ICON));
    popup.addMenuItem(new ModernTwoLineMenuItem("TSS", "Annotate for gene TSS within a region", ICON));

    // popup.addMenuItem(new ModernMenuSeparator());

    popup.addMenuItem(
        new ModernMenuHelpItem("Help with annotating regions...", "geneannotation.help.url").setTextOffset(48));

    mButton = new RibbonLargeDropDownButton(ICON, popup);
    mButton.setToolTip("Annotate", "Annotate genomic regions.");

    mWindow.getRibbon().getToolbar("Genomic").getSection("Annotation").add(mButton);

    mButton.addClickListener(this);
  }

  @Override
  public void clicked(ModernClickEvent e) {
    if (e.getMessage().equals("Genes")) {
      try {
        annotate();
      } catch (IOException e1) {
        e1.printStackTrace();
      }
    } else if (e.getMessage().equals("TSS")) {
      try {
        tssAnnotate();
      } catch (IOException e1) {
        e1.printStackTrace();
      }
    } else {
      // Do nothing
    }
  }

  private void annotate() throws IOException {
    if (mWindow.getCurrentMatrix() == null) {
      return;
    }

    AnnotationDialog dialog = new AnnotationDialog(mWindow);

    dialog.setVisible(true);

    if (dialog.getStatus() == ModernDialogStatus.CANCEL) {
      return;
    }

    List<String> genomes = dialog.getGenomes();

    if (genomes.size() == 0) {
      ModernMessageDialog.createWarningDialog(mWindow, "You must select at least one annotation set.");

      annotate();

      return;
    }

    AnnotateTask task = new AnnotateTask(mWindow, genomes, dialog.getShowOverlappingGenes(), dialog.getClosestList(),
        dialog.getExt5p(), dialog.getExt3p());

    // task.execute();

    task.doInBackground();
    task.done();
  }

  private void tssAnnotate() throws IOException {
    if (mWindow.getCurrentMatrix() == null) {
      return;
    }

    TSSAnnotationDialog dialog = new TSSAnnotationDialog(mWindow);

    dialog.setVisible(true);

    if (dialog.getStatus() == ModernDialogStatus.CANCEL) {
      return;
    }

    List<String> genomes = dialog.getGenomes();

    TSSAnnotateTask task = new TSSAnnotateTask(mWindow, genomes, dialog.getExt5p(), dialog.getExt3p());

    // task.execute();

    task.doInBackground();
    task.done();
  }
}
