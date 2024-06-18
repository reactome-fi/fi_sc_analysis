- The Mouse_cell_markers.txt was downloaded from http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html on February 15, 2022.

- MouseReactomePathways_Hierarchy_Rel_79_122921.xml: Created by method testPathwayHierarchy() in org.reactome.r3.fi.FINetworkResourceTests in the fiviz_corews project. The XML was copied into a file and then saved manually. 
HumanReactomePathways_Hierachy_Rel_79.xml: Created similiarily as the mouse hierarchy file. The database used a production server download file called release_current_010722.sql at the SSD drive. This file was created on Sept 1, 2023.  

- h.all.v2023.1.Hs.symbols.gmt: Downloaded from MsiDB via https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2023.1.Hs/h.all.v2023.1.Hs.symbols.gmt.

- Hg38-blacklist.v2.bed: used by scenic+ via pyCisTopic. Downloaded from pyCisTopic's GitHub: https://github.com/aertslab/pycisTopic/tree/master/blacklist.

- Pathway_List_In_Hirerarchy_Release_79.txt: This file was produced by running a test method, testGetHierarchicalOrderedPathways(), in org.reactome.immport.ws.test in project immport-ws. The output was in the console and manually saved into this file. Note: Need to check how pathways are arranged in this file. The top level pathways are ordered alphabetically. Other pathways are ordered in width first + depth first. This may need to be changed!

- HOM_MouseHumanSequence.rpt.txt: downloaded from http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt. This file is used to map mouse genes to human genes for CellPhoneDB.

- human_protein-coding_gene.txt is downloaded from https://www.genenames.org/download/statistics-and-files/ for all human potein coding genes.

**Note**: The files in this folder are copied from a private GitHub repo, single-cell-analysis.
