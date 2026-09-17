# Exporting the NIST Mass Spectral Library to MSP Format

NIST EI spectral libraries can be exported to MSP format using Lib2NIST or NIST MS Search, depending on the source-library format.

## Method 1: Convert an Agilent/HP `.L` library with Lib2NIST

1. Open **Lib2NIST** and select **Add Input Libraries/Files** to add the `.L` library.
2. Select the library in the input list, choose **NIST text format (`.msp`)** as the output format, and specify the output folder.
3. In **Options**, set the m/z multiplication factor to **1** and the term added before rounding to **0**.
4. Click **Convert** to generate the MSP file.

## Method 2: Export in batches from NIST MS Search

1. Open **NIST MS Search**, select **Other Search → ID Number**, and choose **mainlib**.
2. Enter a consecutive range of library IDs and retrieve the spectra for that batch.
3. Select all spectra in the batch, right-click, and choose **Send To → Spec List**.
4. In the **Librarian** tab, select the batch spectra and click **Export**. Choose **NIST Text (`.msp`)** and save the file as, for example, `mainlib_part01.MSP`.
5. Repeat for the remaining ID ranges, using a separate file for each batch. Export only the current batch each time, covering the library without gaps or overlap.
6. Select **replib** and repeat the same steps to export the replicate library.
