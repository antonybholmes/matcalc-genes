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
import org.jebtk.bioinformatics.genomic.GenomicEntity;
import org.jebtk.bioinformatics.genomic.GenomicRegion;
import org.jebtk.bioinformatics.genomic.GenomicType;
import org.jebtk.bioinformatics.genomic.Strand;
import org.jebtk.core.Mathematics;
import org.jebtk.core.collections.CollectionUtils;
import org.jebtk.core.collections.UniqueArrayList;
import org.jebtk.core.io.Io;
import org.jebtk.core.sys.SysUtils;
import org.jebtk.core.text.Join;
import org.jebtk.core.text.TextUtils;
import org.jebtk.math.matrix.DataFrame;

import edu.columbia.rdf.matcalc.MainMatCalcWindow;
import edu.columbia.rdf.matcalc.bio.AnnotationService;

public class AnnotateTask extends SwingWorker<Void, Void> {

  private DataFrame mNewModel;
  // private boolean mShow2ndClosest;
  // private boolean mShow1stClosest;
  private boolean mShowOverlappingGenes;
  // private boolean mShow3rdClosest;
  // private boolean mShow4thClosest;
  // private boolean mShow5thClosest;

  private int mExt5p = 5000;
  private int mExt3p = 4000;
  private List<Integer> mClosestGenes;
  private MainMatCalcWindow mWindow;
  private Genome mGenome;

  private static final String INTERGENIC = "intergenic";
  private static final String EXONIC = "exonic";
  private static final String INTRONIC = "intronic";
  private static final String PROMOTER = "promoter";

  public AnnotateTask(MainMatCalcWindow window, Genome genome, boolean allSearch, List<Integer> closestGenes, int ext5p,
      int ext3p) {
    mWindow = window;
    mGenome = genome;
    mShowOverlappingGenes = allSearch;
    mClosestGenes = closestGenes;

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

    // See how many extra columns are needed
    int extra = 0;

    List<String> geneIdTypes = AnnotationService.getInstance().getGeneIdTypes(mGenome, GenomicType.TRANSCRIPT);

    int closestSize = 4 + geneIdTypes.size();

    if (mShowOverlappingGenes) {
      extra += 4 + geneIdTypes.size();
    }

    for (int i = 0; i < mClosestGenes.size(); ++i) {
      extra += closestSize;
    }

    DataFrame matrix = DataFrame.createDataFrame(model.getRows(), model.getCols() + extra);

    DataFrame.copyColumnHeaders(model, matrix);

    //
    // Create the new header
    //

    int c = model.getCols();

    geneIdTypes = AnnotationService.getInstance().getGeneIdTypes(mGenome, GenomicType.TRANSCRIPT);

    if (mShowOverlappingGenes) {
      for (String name : geneIdTypes) {
        matrix.setColumnName(c++, "Overlap " + name);
      }

      matrix.setColumnName(c++, "Strand");

      matrix.setColumnName(c++, "Region Labels");

      matrix.setColumnName(c++, "Region Relative To Gene (prom=-" + (mExt5p / 1000) + "/+" + (mExt3p / 1000) + "kb)");

      matrix.setColumnName(c++, "Region TSS Distance (bp)");
    }

    for (int closest : mClosestGenes) {
      c = createClosestHeader(TextUtils.formatOrder(closest), geneIdTypes, c, matrix);
    }

    //
    // Process each row
    //

    for (int i = 0; i < model.getRows(); ++i) {
      matrix.copyRow(model, i, i);

      c = model.getCols();

      GenomicRegion region = null;

      if (Io.isEmptyLine(model.getText(i, 0))) {
        region = null;
      } else if (model.getText(i, 0).contains(TextUtils.NA)) {
        region = null;
      } else if (GenomicRegion.isGenomicRegion(model.getText(i, 0))) {
        region = GenomicRegion.parse(mGenome, model.getText(i, 0));
      } else {
        // three column format

        region = new GenomicRegion(ChromosomeService.getInstance().chr(mGenome, model.getText(i, 0)),
            TextUtils.parseInt(model.getText(i, 1)), TextUtils.parseInt(model.getText(i, 2)));
      }

      // Skip empty annotation
      if (region == null) {
        continue;
      }

      GenomicRegion midRegion = GenomicRegion.midRegion(region);

      // start == end
      int midPoint = midRegion.getStart();

      //
      // If we want overlapping genes
      //

      if (mShowOverlappingGenes) {
        geneIdTypes = AnnotationService.getInstance().getGeneIdTypes(mGenome, GenomicType.TRANSCRIPT);

        GenesDB tssSearch = AnnotationService.getInstance().getSearch(mGenome, mExt5p, mExt3p);

        List<GenomicElement> overlappingResults = tssSearch.overlapping(mGenome, midRegion, GenomicType.TRANSCRIPT, 1);

        // Test midpoint
        // List<GenomicElement> overlappingResults = tssSearch.overlapping(mGenome,
        // midRegion,
        // GenomicType.TRANSCRIPT);

        int n = overlappingResults.size();

        // List<String> allRefseqs = new ArrayList<String>();
        Set<String> labels = new HashSet<String>(n);
        List<String> allClassifications = new ArrayList<String>(n);
        List<String> allTssDist = new ArrayList<String>(n);
        // List<String> allGenes = new ArrayList<String>(n);
        List<Character> allStrands = new ArrayList<Character>(n);

        for (GenomicElement gene : overlappingResults) {
          // Select one refseq
          // String refseq = CollectionUtils.sort(entrezMap.get(gene)).get(0);

          // allGenes.add(symbolMap.get(refseq).getSymbol());
          allStrands.add(Strand.toChar(gene.getStrand()));

          // allRefseqs.add(refseq);

          Set<String> classifications = new HashSet<String>();

          classify((GenomicEntity) gene, region, midPoint, classifications, allTssDist);

          String label = Join.onComma().values(CollectionUtils.reverse(CollectionUtils.sort(classifications)))
              .toString();

          allClassifications.add(label);

          // Add all the unique labels we find (simplified classification).
          labels.add(label);
        }

        if (overlappingResults.size() > 0) {
          for (String name : geneIdTypes) {
            matrix.set(i, c++, formatGeneIds(name, overlappingResults));
          }

          matrix.set(i, c++, TextUtils.scJoin(allStrands));
          matrix.set(i, c++, TextUtils.scJoin(CollectionUtils.sort(labels)));
          matrix.set(i, c++, TextUtils.scJoin(allClassifications));
          matrix.set(i, c++, TextUtils.scJoin(allTssDist));
        } else {
          // Fill the blanks
          c = repeatNA(geneIdTypes.size() + 1, i, c, matrix);
          matrix.set(i, c++, INTERGENIC);
          matrix.set(i, c++, INTERGENIC);
          matrix.set(i, c++, TextUtils.NA);
        }
      }

      //
      // Handle the closest genes
      //

      if (mClosestGenes.size() > 0) {

        // We select the maximum nth closest gene desired
        int maxClosest = Mathematics.maxInt(mClosestGenes); // - 1;

        geneIdTypes = AnnotationService.getInstance().getGeneIdTypes(mGenome, GenomicType.TRANSCRIPT);

        // BinarySearch<AnnotationGene> tssSearch = AnnotationService
        // .getInstance().getBinarySearch(mGenome, mExt5p, mExt3p);

        GenesDB tssSearch = AnnotationService.getInstance().getSearch(mGenome, mExt5p, mExt3p);

        // System.err.println("closest " + closest + " " + c + " " +
        // tssSearch.getName());

        // Find all closest genes between 1 and n
        List<List<GenomicElement>> results = tssSearch.nthClosest(mGenome, midRegion, maxClosest,
            GenomicType.TRANSCRIPT);

        for (int closest : mClosestGenes) {
          // System.err.println(region.getLocation() + ":" + gene.getSymbol()
          // + " " + gene.getLocation());

          // Select the ones of interest
          c = annotateClosest(results.get(closest - 1), geneIdTypes, region, midPoint, i, c, matrix);
        }
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

  private int annotateClosest(final List<GenomicElement> results, final List<String> geneIdTypes, GenomicRegion region,
      int midPoint, int row, int c, DataFrame matrix) {

    /*
     * Map<String, AnnotationGene> entrezMap = new HashMap<String,
     * AnnotationGene>();
     * 
     * if (results != null) { for (AnnotationGene gene : results) {
     * //System.err.println(gene.getRefSeq() + " " + gene.getSymbol());
     * 
     * entrezMap.put(gene.getRefSeq(), gene); } }
     */

    if (results.size() > 0) {
      // If there are more than one, pick the first
      GenomicElement gene = results.get(0); // CollectionUtils.sortKeys(entrezMap).get(0);

      Set<String> classifications = new HashSet<String>();

      int tssDist = classify((GenomicEntity) gene, region, midPoint, classifications, new ArrayList<String>());

      /*
       * int tssDist = Integer.MAX_VALUE;
       * 
       * //for (AnnotationGene gene : results) { tssDist = Math.min(tssDist,
       * AnnotationGene.getTssMidDist(gene, midPoint));
       * 
       * if (withinBounds(gene, midPoint)) { if (tssDist >= -mExt5p && tssDist <=
       * mExt3p) { classifications.add(PROMOTER); }
       * 
       * boolean inExon = false;
       * 
       * for (GenomicRegion exon : gene.getExons()) { if
       * (GenomicRegion.within(midPoint, exon)) { classifications.add(EXONIC); inExon
       * = true; break; } }
       * 
       * if (!inExon && midPoint >= gene.getRegion().getStart() && midPoint <=
       * gene.getRegion().getEnd()) { classifications.add(INTRONIC); } } //}
       * 
       * if (classifications.contains(EXONIC) && classifications.contains(INTRONIC)) {
       * classifications.remove(EXONIC); }
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

  private int classify(GenomicEntity gene, GenomicRegion region, int midPoint, Set<String> classifications,
      List<String> allTssDist) {
    int tssDist = GenomicElement.getTssMidDist(gene, midPoint);

    if (withinBounds(gene, midPoint)) {
      if (tssDist >= -mExt5p && tssDist <= mExt3p) {
        classifications.add(PROMOTER);
      }

      boolean inExon = false;

      for (GenomicRegion exon : gene.getChildren(GenomicType.EXON)) {
        if (GenomicRegion.within(midPoint, exon)) {
          classifications.add(EXONIC);
          inExon = true;
          break;
        }
      }

      if (!inExon && midPoint >= gene.getStart() && midPoint <= gene.getEnd()) {
        classifications.add(INTRONIC);
      }
    }
    // }

    if (classifications.contains(EXONIC) && classifications.contains(INTRONIC)) {
      classifications.remove(EXONIC);
    }

    if (classifications.size() == 0) {
      classifications.add(INTERGENIC);
    }

    if (gene.getSymbol().equals("MYBPHL")) {
      SysUtils.err().println("aha", tssDist, withinBounds(gene, midPoint), gene, region, classifications);
    }

    /*
     * if (classifications.contains(PROMOTER)) {
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
