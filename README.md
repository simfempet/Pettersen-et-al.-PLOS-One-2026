# Pettersen-et-al.-PLOS-One-2026
## Distinct inflammatory and barrier-disruptive podocyte injury programs revealed by transcriptomics and glomerulus-on-chip modeling

### Abstract
Podocytes are central in the organization of the glomerular filtration barrier and function as a major signaling hub for cell-to-cell communication within the glomerulus. Disruption of podocyte biology is therefore central to the progression of chronic kidney disease. To this end, in vitro systems that faithfully resemble healthy podocyte biology and disease perturbation are needed to dissect mechanisms driving podocyte injury and to enable translational studies. Therefore, we aimed to perform a detailed characterization of conditionally immortalized human podocytes during differentiation and to induce injury by lipopolysaccharide (LPS) or the Src-family kinase inhibitor  PP2. To this end, we combined bulk RNA sequencing, quantitative imaging, and a dynamic glomerulus-on-chip model comprising podocytes and endothelial cells. We observed time-dependent transcriptional activation of podocyte marker genes and increased phosphorylation of  nephrin as key differentiation events. LPS treatment increased the production of pro-inflammatory cytokines, and this closely mirrored immune-related processes observed in patient datasets for diabetic-, hypertensive-, and IgA nephropathy, as well as minimal change disease. PP2 treatment resulted in a marked decrease in nephrin phosphorylation, disruption of transcriptional programs controlling slit diaphragm integrity, and increased barrier permeability assessed using glomerulus-on-chip. Taken together, this study presents a versatile tool for mechanistic studies of glomerular diseases, setting the stage for evaluation of therapies targeting glomerular inflammation and barrier dysfunction.

### Listed files
- Differentiation_DESeq2.R
  - Differential expression comparing undifferentiated cells to fully differentiated cells
  - Gene Ontology enrichment analysis of DEGs
  - Fuzzy c-means clustering of genes along the differentiation trajectory
  - Transcription factor activity inference

 - Treatments_DESeq2.R
   - Differential expression comparing LPS-treated cells vs. controls, and PP2-treated cells vs. controls
   - Gene Ontology enrichment analysis of DEGs
   - Characterization of toxicity response using GSEA
   - Podocyte barrier score calculation
  
 - Proteomic_validation.R
   - Processing of MaxQuant output
   - Data quality control
   - Validation of transcriptome findings

 - Semantic_similarity.R
   - Differential expression analysis of human glomerular disease datasets
   - Gene Ontology enrichment analysis of DEGs
   - Semantic similarity calculation between human glomerular diseases and in vitro results

 - PCA.R
   - Principal component analysis of differentiation and treatment-control samples

- Cell_viability.R
   - Analysis of resazurin sodium salt viability assay

- Nuclei_morphometrics.R
   - Analysis of cell death induction using nuclei morphology
 
- P_NPHS1.R
   - Quantitative analysis of NPHS1-pY1176/1193 immunofluorescence images
 
- Permeability.R
   - Analysis of glomerulus-on-chip permeability assay
   
   
   
