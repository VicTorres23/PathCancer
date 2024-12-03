from flask import Flask, request, render_template, redirect, url_for, send_from_directory
import os
from werkzeug.utils import secure_filename
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import time
import zipfile
import vcf
import networkx as nx
import subprocess

app = Flask(__name__)

UploadFolder = 'uploads/'
DownloadsFolder = 'static/downloads/'
ImagesFolder = 'static/images/'
LogoFolder = 'static/logo/'
TSVFolder = 'TSV_files/'
AllowedExtensions = {'maf', 'vcf'}
app.config['UploadFolder'] = UploadFolder

os.makedirs(UploadFolder, exist_ok=True)
os.makedirs(DownloadsFolder, exist_ok=True)
os.makedirs(ImagesFolder, exist_ok=True)
os.makedirs(LogoFolder, exist_ok=True)
os.makedirs(TSVFolder, exist_ok=True)

LogoPath = os.path.join(LogoFolder, "PathCancer_logo.png")

def VerifyFile(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in AllowedExtensions

def GetCenter(vcf_filename):
    with open(vcf_filename) as vcf_file:
        for line in vcf_file:
            if line.startswith("##center="):
                return line.strip().split('=')[1].strip('"')
    return "Unknown"

def GetNCBIBuild(vcf_filename):
    with open(vcf_filename) as vcf_file:
        for line in vcf_file:
            if line.startswith("##reference="):
                reference = line.strip().split('=')[1]
                if "GRCh37" in reference or "37" in reference:
                    return "37"
                elif "GRCh38" in reference or "38" in reference:
                    return "38"
    return "Unknown"

def GetTumorSampleBarcode(vcf_filename):
    with open(vcf_filename) as vcf_file:
        for line in vcf_file:
            if line.startswith("##SAMPLE=<ID=TUMOR"):
                SampleInfo = line.strip().split("NAME=")[1]
                barcode = SampleInfo.split(",")[0].split("-")
                return "-".join(barcode[:3])
    return "Unknown"

ConsequenceToClassification = {
    "splice_acceptor_variant": "Splice_Site",
    "splice_donor_variant": "Splice_Site",
    "missense_variant": "Missense_Mutation",
    "synonymous_variant": "Silent",
    "stop_gained": "Nonsense_Mutation",
    "frameshift_variant": "Frame_Shift_Ins",
    "inframe_insertion": "In_Frame_Ins",
    "inframe_deletion": "In_Frame_Del",
    "start_lost": "Translation_Start_Site",
    "stop_lost": "Nonstop_Mutation",
    "3_prime_UTR_variant": "3'UTR",
    "5_prime_UTR_variant": "5'UTR",
    "intron_variant": "Intron",
    "upstream_gene_variant": "5'Flank",
    "downstream_gene_variant": "3'Flank",
    "intergenic_variant": "IGR",
    "non_coding_transcript_exon_variant": "RNA",
    "regulatory_region_variant": "Regulatory"
}

def VCFtoMAF(vcf_filename):
    vcf_file = vcf.Reader(filename=vcf_filename)
    tumor_sample_barcode = GetTumorSampleBarcode(vcf_filename)
    HUGO_SYMBOL = []
    ENTREZ_ID = []
    CHROMOSOME = []
    START_POS = []
    END_POS = []
    variant_classifications = []
    variant_types = []
    reference_alleles = []
    tumor_seq_allele2 = []
    PROTEIN_CHANGE = []
    i_TumorVAF_WU = []
    i_transcript_name = []

    for record in vcf_file:
        fields = record.INFO["CSQ"][0].split("|")
        HUGO_SYMBOL.append(fields[3])
        ENTREZ_ID.append(fields[-2])
        CHROMOSOME.append(record.CHROM.replace("chr", ""))
        START_POS.append(record.POS)
        END_POS.append(record.POS)

        consequences = fields[1].split("&")
        classification = "Unknown"
        for consequence in consequences:
            if consequence in ConsequenceToClassification:
                classification = ConsequenceToClassification[consequence]
                break
        variant_classifications.append(classification)

        ref = record.REF
        alt = record.ALT[0]
        if len(ref) == len(alt) == 1:
            variant_type = "SNP"
        elif len(ref) == 1 and len(alt) > 1:
            variant_type = "INS"
        elif len(ref) > 1 and len(alt) == 1:
            variant_type = "DEL"
        else:
            variant_type = "Complex"
        variant_types.append(variant_type)

        reference_alleles.append(record.REF)
        tumor_seq_allele2.append(record.ALT[0] if record.ALT else None)
        PROTEIN_CHANGE.append(fields[11].split("(")[1].strip(")") if "(" in fields[11] else "N/A")

        tumor_info = record.genotype("TUMOR")["AD"]
        if tumor_info and len(tumor_info) == 2:
            ref_depth, alt_depth = tumor_info
            vaf = (alt_depth / (ref_depth + alt_depth)) * 100 if (ref_depth + alt_depth) > 0 else 0
            i_TumorVAF_WU.append(vaf)
        else:
            i_TumorVAF_WU.append(None)
        i_transcript_name.append(fields[33])

    maf = pd.DataFrame({
        "Hugo_Symbol": HUGO_SYMBOL,
        "Entrez_Gene_Id": ENTREZ_ID,
        "Center": GetCenter(vcf_filename),
        "NCBI_Build": GetNCBIBuild(vcf_filename),
        "Chromosome": CHROMOSOME,
        "Start_Position": START_POS,
        "End_Position": END_POS,
        "Strand": ["+"] * len(HUGO_SYMBOL),
        "Variant_Classification": variant_classifications,
        "Variant_Type": variant_types,
        "Reference_Allele": reference_alleles,
        "Tumor_Seq_Allele1": reference_alleles,
        "Tumor_Seq_Allele2": tumor_seq_allele2,
        "Tumor_Sample_Barcode": [tumor_sample_barcode] * len(HUGO_SYMBOL),
        "Protein_Change": PROTEIN_CHANGE,
        "i_TumorVAF_WU": i_TumorVAF_WU,
        "i_transcript_name": i_transcript_name
    })
    maf = maf[maf['Hugo_Symbol'] != '']
    return maf

def getZipFile(HTMLpath, imagespath, LogoPath, ZIPpath, data_files):
    with zipfile.ZipFile(ZIPpath, "w") as zip:
        zip.write(HTMLpath, os.path.basename(HTMLpath))
        for image_path in imagespath:
            zip.write(image_path,
                      os.path.join('images', os.path.basename(image_path)))
        zip.write(LogoPath, os.path.join('logo', os.path.basename(LogoPath)))
        for data_file in data_files:
            zip.write(data_file, os.path.join('data', os.path.basename(data_file)))

def ExtractMafStatistics(df):

    maf = df

    NumSamples = maf['Tumor_Sample_Barcode'].nunique() if 'Tumor_Sample_Barcode' in maf else 0
    NumVariants = len(maf)

    ReferenceGenome = maf['NCBI_Build'].iloc[0] if 'NCBI_Build' in maf else "Unknown"

    AnnotationResources = maf['Center'].iloc[0] if 'Center' in maf else "Unknown"

    MissenseCount = maf[maf['Variant_Classification'] == 'Missense_Mutation'].shape[0]
    NonsenseCount = maf[maf['Variant_Classification'] == 'Nonsense_Mutation'].shape[0]
    SilentCount = maf[maf['Variant_Classification'] == 'Silent'].shape[0]

    statistics = {
        'Number of Samples': NumSamples,
        'Number of Variants': NumVariants,
        'Reference Genome': ReferenceGenome,
        'Annotation Resources': AnnotationResources,
        'Missense Variants': MissenseCount,
        'Nonsense Variants': NonsenseCount,
        'Silent Variants': SilentCount
    }

    return statistics

def getClassificationPlot(df, OutputPath):

    maf = df

    sns.set(style="whitegrid")

    plt.figure(figsize=(8, 6))
    variant_classification_counts = maf['Variant_Classification'].value_counts()
    sns.barplot(x=variant_classification_counts.values, y=variant_classification_counts.index, palette="Blues_r")
    plt.title("Variant Classification")
    plt.xlabel("Count")
    plt.ylabel("Classification")
    plt.tight_layout()
    output_path = os.path.join(OutputPath, "variant_classification.png")
    plt.savefig(output_path)
    plt.close()

    return output_path

def getMutationNetwork(CoOcurrenceMatrix, output_path):

    G = nx.Graph()

    for gene1 in CoOcurrenceMatrix.index:
        for gene2 in CoOcurrenceMatrix.columns:
            if gene1 != gene2 and CoOcurrenceMatrix.at[gene1, gene2] > 0:
                G.add_edge(gene1, gene2, weight=CoOcurrenceMatrix.at[gene1, gene2])

    edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())

    plt.figure(figsize=(10, 10))
    pos = nx.spring_layout(G, seed=42)
    nx.draw_networkx_nodes(G, pos, node_size=700, node_color='lightblue')
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight="bold")
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=np.log2(weights) + 1, alpha=0.6)
    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels={(gene1, gene2): str(weight) for gene1, gene2, weight in G.edges(data='weight')}
    )

    plt.title("Network of Mutation Relationships", fontsize=15)
    plt.axis('off')

    plt.savefig(output_path, dpi=300)
    plt.close()

def getVariantTypePlot(df, OutputPath):

    maf = df

    sns.set(style="whitegrid")
    plt.figure(figsize=(8, 6))
    variant_type_counts = maf['Variant_Type'].value_counts()

    sns.barplot(
        x=variant_type_counts.values,
        y=variant_type_counts.index,
        palette="Reds_r"
    )

    plt.title("Variant Type")
    plt.xlabel("Count")
    plt.ylabel("Type")
    plt.tight_layout()

    output_path = os.path.join(OutputPath, "variant_type.png")

    plt.savefig(output_path)
    plt.close()

    return output_path

def GenerateGraphs(FilePath):

    _, FileExtension = os.path.splitext(FilePath)
    FileExtension = FileExtension.lower()

    OutputPath = ImagesFolder
    DataFolder = TSVFolder

    if FileExtension == '.maf':
        df = pd.read_csv(FilePath, delimiter='\t')
        Index = "Tumor_Sample_Barcode"
        Columns = "Hugo_Symbol"
    elif FileExtension == '.vcf':
        VCFFiles = [os.path.join(UploadFolder, f) for f in os.listdir(UploadFolder) if f.endswith('.vcf')]
        MAFdfs = [VCFtoMAF(vcf_file) for vcf_file in VCFFiles]
        df = pd.concat(MAFdfs, ignore_index=True)
        Index = "Tumor_Sample_Barcode"
        Columns = "Hugo_Symbol"
    else:
        return 'Unsupported file type', None

    pivot_table = pd.pivot_table(df, index=Index, columns=Columns, aggfunc='size', fill_value=0)

    TopGenes = pivot_table.sum(axis=0).nlargest(20).index
    FilteredPivotTable = pivot_table[TopGenes]

    CoOcurrenceMatrix = np.dot(FilteredPivotTable.T, FilteredPivotTable)
    CoOcurrenceMatrix = pd.DataFrame(
        CoOcurrenceMatrix,
        index=FilteredPivotTable.columns,
        columns=FilteredPivotTable.columns
    )

    PivotTablePath = os.path.join(DataFolder, "pivot_table.tsv")
    pivot_table.to_csv(PivotTablePath, sep='\t')

    CoOcurrenceTSV = os.path.join(DataFolder, "co_occurrence_matrix.tsv")
    CoOcurrenceMatrix.to_csv(CoOcurrenceTSV, sep='\t')

    CoOcurrencePath = os.path.join(ImagesFolder, 'co_occurrence_heatmap.png')
    plt.figure(figsize=(6, 5))
    sns.heatmap(CoOcurrenceMatrix, annot=True, cmap="Blues",
                cbar_kws={'label': 'Co-Occurrence Count'}, square=True, xticklabels=True, yticklabels=True)
    plt.title("Top 20 Gene Mutation Co-Occurrence Heatmap")
    plt.xlabel("Gene")
    plt.ylabel("Gene")
    plt.savefig(CoOcurrencePath)
    plt.close()

    MutationNetworkPath = os.path.join(ImagesFolder, 'mutation_network.png')
    getMutationNetwork(CoOcurrenceMatrix, MutationNetworkPath)

    MutualExclusivityPath = os.path.join(ImagesFolder, 'mutual_exclusivity_heatmap.png')
    MutualExclusivityMatrix = pd.DataFrame(0.0, index=TopGenes, columns=TopGenes)
    for gene1 in TopGenes:
        for gene2 in TopGenes:
            if gene1 != gene2:
                ExclusiveCount = np.sum((FilteredPivotTable[gene1] + FilteredPivotTable[gene2]) == 1)
                MutualExclusivityMatrix.loc[gene1, gene2] = ExclusiveCount / len(FilteredPivotTable)

    MutualExclusivityTSV = os.path.join(DataFolder, "mutual_exclusivity_matrix.tsv")
    MutualExclusivityMatrix.to_csv(MutualExclusivityTSV, sep='\t')

    VariantDataPath = os.path.join(DataFolder, "variant_data.tsv")
    df.to_csv(VariantDataPath, sep='\t', index=False)

    plt.figure(figsize=(6, 5))
    sns.heatmap(MutualExclusivityMatrix, annot=True, fmt=".1f", cmap="Reds", cbar_kws={'label': 'Mutual Exclusivity Score'},
                square=True, xticklabels=True, yticklabels=True, annot_kws={"size": 8})
    plt.title("Top 20 Gene Mutation Mutual Exclusivity Heatmap")
    plt.xlabel("Gene")
    plt.ylabel("Gene")
    plt.savefig(MutualExclusivityPath)
    plt.close()

    VariantClassificationPath = getClassificationPlot(df, OutputPath)
    VariantTypePath = getVariantTypePlot(df, OutputPath)

    ScriptDir = os.path.dirname(os.path.abspath(__file__))
    rScriptPath = os.path.join(ScriptDir, "Pathway_Analysis.R")
    subprocess.run(["C:/Program Files/R/R-4.3.3/bin/Rscript", rScriptPath, FilePath, OutputPath], check=True)

    PathwayBarplotPath = os.path.join(OutputPath, "pathway_barplot.png")

    MAFStats = ExtractMafStatistics(df)

    return (
        CoOcurrencePath,
        MutualExclusivityPath,
        MutationNetworkPath,
        [PathwayBarplotPath, VariantClassificationPath], VariantTypePath, MAFStats,
        [PivotTablePath, CoOcurrencePath, MutualExclusivityPath, VariantDataPath]
    )

@app.route('/', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        if 'file' not in request.files:
            return 'No file part in the request', 400
        file = request.files['file']
        if file.filename == '':
            return 'No selected file', 400
        if file and VerifyFile(file.filename):

            filename = secure_filename(file.filename)
            FilePath = os.path.join(app.config['UploadFolder'], filename)
            file.save(FilePath)

            (CoOcurrencePath, MutualExclusivityPath, MutationNetworkPath, r_generated_images, VariantTypePath, MAFStats, data_files) = GenerateGraphs(FilePath)

            timestamp = int(time.time())

            DownloadHTML = f"pathcancer_report_{timestamp}.html"
            DownloadHTMLPath = os.path.join(DownloadsFolder, DownloadHTML)
            with open(DownloadHTMLPath, "w") as file:
                file.write(render_template(
                    "result.html",
                    co_occurrence_image="images/" + os.path.basename(CoOcurrencePath),
                    mutual_exclusivity_image="images/" + os.path.basename(MutualExclusivityPath),
                    logo_image="logo/PathCancer_logo.png",
                    MAFStats=MAFStats,
                    mutation_network_image="images/" + os.path.basename(MutationNetworkPath),
                    pathway_barplot="images/" + os.path.basename(r_generated_images[0]),
                    variant_classification_image="images/variant_classification.png",
                    variant_type_image="images/" + os.path.basename(VariantTypePath),
                    download_link=True
                ))

            ZIPName = f"pathcancer_report_{timestamp}.zip"
            ZIPPath = os.path.join(DownloadsFolder, ZIPName)
            getZipFile(DownloadHTMLPath, [CoOcurrencePath, MutualExclusivityPath, MutationNetworkPath] + r_generated_images, LogoPath,
                       ZIPPath, data_files)

            return render_template('ready.html', download_link=True, download_zip=ZIPName)
    return render_template('index.html')

@app.route("/download/<filename>")
def download_file(filename):
    return send_from_directory(DownloadsFolder, filename, as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)