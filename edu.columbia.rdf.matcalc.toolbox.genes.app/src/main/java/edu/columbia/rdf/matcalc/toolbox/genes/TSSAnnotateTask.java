package edu.columbia.rdf.matcalc.toolbox.genes;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.SwingWorker;

import org.jebtk.bioinformatics.genomic.ChromosomeService;
import org.jebtk.bioinformatics.genomic.GenesDB;
import org.jebtk.bioinformatics.genomic.Genome;
import org.jebtk.bioinformatics.genomic.GenomicElement;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.GenomicType;
import org.jebtk.bioinformatics.genomic.Strand;
import org.jebtk.core.collections.ArrayListCreator;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.DefaultTreeMap;
import org.jebtk.core.collections.DefaultTreeMapCreator;
import org.jebtk.core.collections.IterMap;
import org.jebtk.core.collections.UniqueArrayList;
import org.jebtk.core.io.Io;
import org.jebtk.core.text.Join;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.matrix.DataFrame;

import edu.columbia.rdf.matcalc.MainMatCalcWindow;
import edu.columbia.rdf.matcalc.bio.AnnotationService;

public class TSSAnnotateTask extends SwingWorker<Void, Void> {

  private DataFrame mNewModel;
  // private boolean mShow2ndClosest;
  // private boolean mShow1stClosest;
  // private boolean mShow3rdClosest;
  // private boolean mShow4thClosest;
  // private boolean mShow5thClosest;

  private int mExt5p = 5000;
  private int mExt3p = 4000;
  private MainMatCalcWindow mWindow;
  private Genome mGenome;

  public TSSAnnotateTask(MainMatCalcWindow window, Genome genome, int ext5p, int ext3p) {
    mWindow = window;
    mGenome = genome;

    mExt5p = ext5p;
    mExt3p = ext3p;
  }

  @Override
  public Void doInBackground() {
    try {
      mNewModel = annotation();
    } catch (Exception e) {
      e.printStackTrace();
    }

    return null;
  }

  @Override
  public void done() {
    if (mNewModel != null) {
      mWindow.history().addToHistory("Add Annotation", mNewModel);
    }
  }

  private DataFrame annotation() throws IOException, ParseException {
    DataFrame model = mWindow.getCurrentMatrix();

    //
    // Process each row
    //

    IterMap<GenomicRegion, IterMap<Genome, List<GenomicElement>>> overlappingResults = DefaultTreeMap
        .create(new DefaultTreeMapCreator<Genome, List<GenomicElement>>(new ArrayListCreator<GenomicElement>(100)));

    List<GenomicRegion> regions = new ArrayList<GenomicRegion>();

    for (int i = 0; i < model.getRows(); ++i) {
      GenomicRegion region = null;

      if (Io.isEmptyLine(model.getText(i, 0))) {
        region = null;
      } else if (model.getText(i, 0).contains(TextUtils.NA)) {
        region = null;
      } else if (GenomicRegion.isGenomicRegion(model.getText(i, 0))) {
        region = GenomicRegion.parse(Genome.HG19, model.getText(i, 0));
      } else {
        // three column format
        region = new GenomicRegion(ChromosomeService.getInstance().chr(Genome.HG19, model.getText(i, 0)),
            TextUtils.parseInt(model.getText(i, 1)), TextUtils.parseInt(model.getText(i, 2)));
      }

      // Skip empty annotation
      if (region == null) {
        continue;
      }

      GenomicRegion newRegion = GenomicRegion.create(region.mChr, region.mStart + mExt5p, region.mEnd - mExt3p);

      regions.add(newRegion);

      GenesDB tssSearch = AnnotationService.getInstance().getSearch(mGenome);

      List<GenomicElement> results = tssSearch.find(mGenome, newRegion, GenomicType.TRANSCRIPT, 1);

      for (GenomicElement gene : results) {
        // System.err.println(region.getLocation() + ":" + gene.getSymbol() +
        // " " +
        // gene);

        // Check that Gene tss lies within region
        if (GenomicRegion.overlap(newRegion, gene.getTss()) != null) {
          overlappingResults.get(newRegion).get(mGenome).add(gene);
        }
      }
    }

    // How many extra columns we need

    // distance from 5' and 3' end plus strand
    int extra = 0;

    List<String> geneIdTypes = AnnotationService.getInstance().getGeneIdTypes(mGenome, GenomicType.TRANSCRIPT);

    extra += geneIdTypes.size() + 4;

    DataFrame matrix = DataFrame.createDataFrame(model.getRows(), model.getCols() + extra);

    DataFrame.copyColumnHeaders(model, matrix);

    // Add the extra columns
    int c = model.getCols();

    geneIdTypes = AnnotationService.getInstance().getGeneIdTypes(mGenome, GenomicType.TRANSCRIPT);

    for (String name : geneIdTypes) {
      matrix.setColumnName(c++, mGenome + " TSS " + name);
    }

    matrix.setColumnName(c++, mGenome + " strand");
    matrix.setColumnName(c++, mGenome + " TSS");
    matrix.setColumnName(c++, mGenome + " TSS 5' distance");
    matrix.setColumnName(c++, mGenome + " TSS 3' distance");

    // Copy existing rows
    for (int i = 0; i < model.getRows(); ++i) {
      matrix.copyRow(model, i, i);

      GenomicRegion region = regions.get(i);

      c = model.getCols();

      geneIdTypes = AnnotationService.getInstance().getGeneIdTypes(mGenome, GenomicType.TRANSCRIPT);

      List<GenomicElement> results = overlappingResults.get(region).get(mGenome);

      int n = results.size();

      if (n > 0) {
        List<GenomicRegion> allTss = new ArrayList<GenomicRegion>(n);
        List<Integer> dist5p = new ArrayList<Integer>(n);
        List<Integer> dist3p = new ArrayList<Integer>(n);
        List<Character> allStrands = new ArrayList<Character>(n);

        for (GenomicElement gene : results) {
          int pd5;
          int pd3;

          if (Strand.isSense(gene.getStrand())) {
            pd5 = gene.getTss().getStart() - region.getStart();
            pd3 = gene.getTss().getStart() - region.getEnd();
          } else {
            pd5 = region.getStart() - gene.getTss().getStart();
            pd3 = region.getEnd() - gene.getTss().getStart();
          }

          dist5p.add(pd5);
          dist3p.add(pd3);
          allTss.add(gene.getTss());
          allStrands.add(Strand.toChar(gene.getStrand()));
        }

        for (String name : geneIdTypes) {
          matrix.set(i, c++, formatGeneIds(name, results));
        }

        matrix.set(i, c++, TextUtils.scJoin(allStrands));
        matrix.set(i, c++, TextUtils.scJoin(allTss));
        matrix.set(i, c++, TextUtils.scJoin(dist5p));
        matrix.set(i, c++, TextUtils.scJoin(dist3p));
      } else {
        repeatNA(extra, i, c++, matrix);
      }
    }

    return matrix;
  }

  /**
   * Repeat 'n/a' in a matrix cell.
   * 
   * @param n      The number of 'n/a' to write.
   * @param row    The row to write on.
   * @param c      The starting column.
   * @param matrix The matrix to update.
   * @return The new column (equal to c + n).
   */
  public int repeatNA(int n, int row, int c, DataFrame matrix) {
    for (int i = 0; i < n; ++i) {
      matrix.set(row, c++, TextUtils.NA);
    }

    return c;
  }

  /**
   * Create a semi colon separated list of gene ids from a set of genes.
   * 
   * @param name
   * @param genes
   * @return
   */
  private String formatGeneIds(String name, List<GenomicElement> genes) {
    List<String> items = new UniqueArrayList<String>(genes.size());

    for (GenomicElement gene : genes) {
      items.add(gene.getProperty(name));
    }

    return Join.onSemiColon().values(items).toString();
  }

  private int createClosestHeader(String label, List<String> geneIdTypes, int c, DataFrame matrix) {
    // matrix.setColumnName(c++, label + " Closest RefSeq");
    // matrix.setColumnName(c++, label + " Closest Gene Entrez ID");
    // matrix.setColumnName(c++, label + " Closest Gene Symbol");

    for (String name : geneIdTypes) {
      matrix.setColumnName(c++, label + " Closest " + name);
    }

    matrix.setColumnName(c++, label + " Closest Gene Strand");
    matrix.setColumnName(c++, label + " Closest TSS Location");
    matrix.setColumnName(c++,
        label + " Closest Region Relative To Gene (prom=-" + (mExt5p / 1000) + "/+" + (mExt3p / 1000) + "kb)");
    matrix.setColumnName(c++, label + " Closest Region TSS Distance (bp)");

    return c;
  }

  private int annotateClosest(final List<GenomicElement> results, final List<String> geneIdTypes, int midPoint, int row,
      int c, DataFrame matrix) {

    if (results.size() > 0) {
      // If there are more than one, pick the first
      GenomicElement gene = results.get(0); // CollectionUtils.sortKeys(entrezMap).get(0);

      Set<String> classifications = new HashSet<String>();

      int tssDist = classify(gene, midPoint, classifications, new ArrayList<String>());

      /*
       * int tssDist = Integer.MAX_VALUE;
       * 
       * //for (GenomicElement gene : results) { tssDist = Math.min(tssDist,
       * GenomicElement.getTssMidDist(gene, midPoint));
       * 
       * if (withinBounds(gene, midPoint)) { if (tssDist >= -mExt5p && tssDist <=
       * mExt3p) { classifications.add("promoter"); }
       * 
       * boolean inExon = false;
       * 
       * for (GenomicRegion exon : gene.getExons()) { if
       * (GenomicRegion.within(midPoint, exon)) { classifications.add("exonic");
       * inExon = true; break; } }
       * 
       * if (!inExon && midPoint >= gene.getStart() && midPoint <= gene.getEnd()) {
       * classifications.add("intronic"); } } //}
       * 
       * if (classifications.contains("exonic") &&
       * classifications.contains("intronic")) { classifications.remove("exonic"); }
       * 
       * if (classifications.size() == 0) { classifications.add("intergenic"); }
       */

      for (String type : geneIdTypes) {
        matrix.set(row, c++, gene.getProperty(type));
      }

      // matrix.set(row, c++, refseq);
      // matrix.set(row, c++, entrezMap.get(refseq).getEntrez());
      // matrix.set(row, c++, entrezMap.get(refseq).getSymbol());

      matrix.set(row, c++, Strand.toString(gene.getStrand()));
      matrix.set(row, c++, gene.getTss().getLocation());
      matrix.set(row, c++, TextUtils.commaJoin(CollectionUtils.reverse(CollectionUtils.sort(classifications))));
      matrix.set(row, c++, tssDist);
    } else {
      // matrix.set(row, c++, TextUtils.NA);
      // matrix.set(row, c++, TextUtils.NA);
      // matrix.set(row, c++, TextUtils.NA);

      // Fill the blanks
      c = repeatNA(geneIdTypes.size() + 4, row, c, matrix);
    }

    return c;
  }

  private int classify(GenomicElement gene, int midPoint, Set<String> classifications, List<String> allTssDist) {
    int tssDist = GenomicElement.getTssMidDist(gene, midPoint);

    if (withinBounds(gene, midPoint)) {
      if (tssDist >= -mExt5p && tssDist <= mExt3p) {
        classifications.add("promoter");
      }

      boolean inExon = false;

      for (GenomicRegion exon : gene.getChildren(GenomicType.EXON)) {
        if (GenomicRegion.within(midPoint, exon)) {
          classifications.add("exonic");
          inExon = true;
          break;
        }
      }

      if (!inExon && midPoint >= gene.getStart() && midPoint <= gene.getEnd()) {
        classifications.add("intronic");
      }
    }
    // }

    if (classifications.contains("exonic") && classifications.contains("intronic")) {
      classifications.remove("exonic");
    }

    if (classifications.size() == 0) {
      classifications.add("intergenic");
    }

    /*
     * if (classifications.contains("promoter")) {
     * allTssDist.add(Integer.toString(tssDist)); } else {
     * allTssDist.add(TextUtils.NA); }
     */

    // Report distance regardless of whether we are in the boundaries
    // of the promoter region
    allTssDist.add(Integer.toString(tssDist));

    return tssDist;
  }

  /**
   * Returns true if the point is within the gene bounds.
   * 
   * @param gene
   * @param midPoint
   * @return
   */
  private boolean withinBounds(GenomicElement gene, int midPoint) {
    if (gene.getStrand() == Strand.SENSE) {
      return midPoint >= gene.getStart() - mExt5p && midPoint <= gene.getEnd();
    } else {
      return midPoint >= gene.getStart() && midPoint <= gene.getEnd() + mExt5p;
    }
  }
}
