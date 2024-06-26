# AML_Tcell
Single cell RNA-seq with AML and T cell data 

### Step 1: Read and Prepare Data
	 ~ Read the scRNA-seq data files for the AML and T cells 
	 ~ Make sure files are separated into each set of barcode, feature, and matrix files
	 ~ The files are then renamed "barcodes.tsv.gz" , "features.tsv.gz", and "matrix.mtx.gz" to match the format expected by the Read10X function
  	 ~ Merge files into one dataset and rename

### Step 2: Create Seurat Object
        Create Seurat Object for merged dataset. Seurat objects are used to store both the expression data and associated metadata.
     	CreateSeuratObject: Initializes a Seurat object.
              Parameters:
    	        counts: Raw gene expression data.
     	        project: Name for the project.
     	        min.cells: Minimum number of cells expressing a gene for it to be included.
     	        min.features: Minimum number of genes expressed in a cell for it to be included.
              Result: 
                class, features(genes), samples(cells), assay, active assay, layer

### Step 3: Quality Control (QC) 
        Calculate quality control metrics
                Filter out low quality cells and genes and normalize the data
                Calculate percent mitochondrial QC metrics
                QC metrics: "nfeature_RNA", "nCount_RNA", "percent.mt"
                Low quality cells or empty droplets often have very few genes
                Cell doublets or multiplets have high values of nfeature_RNA & nCount_RNA
                Low quality cells often have high percentage of mitochondrial genes
                PercentageFeatureSet: Computes the percentage of mitochondrial gene expression, which is a common QC metric.
<img width="460" alt="3" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/4773ba5e-8b37-4685-b2bd-ccec2b5f4fff">

         Visualize using violin plots
                VlnPlot: Creates violin plots to visualize the distribution of QC metrics.
                
![Picture1](https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/8616bc4f-bd81-42ff-8b44-5e357dd48502)
                

         Filter out low-quality cells
                FeatureScatter: Creates scatter plots to visualize relationships between QC metrics
                ![2](https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/fc733ecd-73b4-490e-adab-27eaa73bab5f)
![2](https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/fc733ecd-73b4-490e-adab-27eaa73bab5f)

### Step 4: Normalization
        Normalize the data to account for differences in sequencing depth across cells.
        NormalizeData: Normalizes the expression data.
        normalization.method = "LogNormalize": Normalizes gene expression values by log-transformation.
        scale.factor = 10000: Scaling factor for normalization.
### Step 5: Identify Highly Variable Features
        Highly variable features are genes whose expression varies significantly across cells, which are often more informative for distinguishing between cell types or conditions.
        FindVariableFeatures: Identifies highly variable genes.
        selection.method = "vst": Variance Stabilizing Transformation method.
        nfeatures = 3000: Number of variable genes to identify.
        VariableFeaturePlot: generates a plot that visualizes the distribution and variability of gene expression across cells.
                look for genes at the right end of x-axis with higher variability in expression across cells
<img width="900" alt="Screenshot 2024-06-25 at 4 47 34 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/5aae8270-658d-4241-b311-2fa6214d92ee">

### Step 6: Scaling Data 
        Centers and scales gene expression values, typically to have zero mean and unit variance.
        ScaleData: Centers and scales the data for each gene.
        features = all.genes: Scales all genes in the dataset.

### Step 7: Perform linear dimensional reduction
        Perform Principal Component Analysis (PCA) on genes
        features = the highly variable genes
        DimPlot and DimHeatmap to look at gene expression and cells along a specific PCA dimension.
<img width="500" alt="Screenshot 2024-06-25 at 4 52 32 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/1852719c-1e0c-4c95-9dae-3ed412d2f853">
<img width="500" alt="Screenshot 2024-06-25 at 4 54 14 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/ea6ace8a-ebdb-47a0-8b3e-a4a310400dd8">
<img width="713" alt="Screenshot 2024-06-25 at 4 56 09 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/a6690ea0-60a2-42d2-a238-76cfe4158368">

### Step 8: Determine Dimensionality 
        Elbow Plot: determining the optimal number of principal components (PC) or dimensions to retain for downstream analysis.
<img width="600" alt="Screenshot 2024-06-25 at 4 57 38 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/d31db91d-5cde-4a26-aa7d-465d81852eb7">

### Step 9: Clustering Cells 
        Cluster the cells by their PCA scores and visualize these clusters using UMAP.
                FindNeighbors: Constructs a k-nearest neighbors graph.
                dims = 1:10: Utilizes the first 10 principal components (PCs).
                FindClusters: Identifies distinct clusters of cells.
                resolution = 0.5: Controls the granularity of the clustering process.
                RunUMAP: Performs Uniform Manifold Approximation and Projection (UMAP) for dimensional reduction and visualization.
                DimPlot: Presents the clustered cells in a 2D representation. 
<img width="712" alt="Screenshot 2024-06-25 at 5 00 24 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/47ac3bd6-bda4-4399-aa12-661b75829fa0">

### Step 10: Cluster biomarkers
        Utilizes FindAllMarkers to identify markers for each cluster against all other cells, focusing on positive markers (upregulated).
<img width="400" alt="Screenshot 2024-06-26 at 1 56 33 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/99485bef-887d-4bf8-9e56-bad3337baff8">
<img width="400" alt="Screenshot 2024-06-26 at 1 58 47 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/414a5993-ed70-48e1-ba3e-3e9ad9269122">
        
        Feature plot: visualize the expression levels of specified features (genes or other variables) across cells in a dataset
<img width="500" alt="Screenshot 2024-06-26 at 2 10 24 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/5cdb024a-6669-4f13-8d11-af36e399ec02">
<img width="500" alt="Screenshot 2024-06-26 at 2 22 26 PM" src="https://github.com/zoey-shawlab/Single-Cell-Seurat/assets/173721950/e43b4232-9375-4aaf-be2d-aa62650e23f7">

        Heatmap: identify and visualize top marker genes (based on differential expression) for each cluster
<img width="900" alt="Screenshot 2024-06-26 at 2 13 33 PM" src="https://github.com/zoey-shawlab/AML_Tcell/assets/173721950/495921de-2f1f-47b5-97dd-2466981af3fd">








        







        



        
