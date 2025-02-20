package cl.util.run;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Iterator;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

public class CombineSolutionCSV {
	public static void main(String[] args) throws IOException
	{
		// Set input folder.
		String inputFolder = "data/solutions";

		// Initialise new StringBuilder.
		int rowIndex = 0;
		try (XSSFWorkbook wbOutput = new XSSFWorkbook()) {
			OutputStream os = new FileOutputStream(new File(inputFolder + "/solution.xlsx"));
			Sheet s = wbOutput.createSheet();
			
			File inputDir = new File(inputFolder);
			for (File f : inputDir.listFiles())
			{
				try (XSSFWorkbook wb = new XSSFWorkbook(f))
				{
					XSSFSheet sheet = wb.getSheetAt(0);  
					Iterator<Row> rowIterator = sheet.rowIterator();
					Row headerRow = rowIterator.next();
					
					if (rowIndex==0)
					{	
						Row newRow = s.createRow(rowIndex++);
						writeCells(newRow, headerRow.cellIterator());
					}
					while (rowIterator.hasNext())
					{
						Row row = rowIterator.next();
						Row newRow = s.createRow(rowIndex++);
						writeCells(newRow, row.cellIterator());
					}
				}
				catch (Exception e)
				{
					// Do nothing.
				}
			}

			// Write results.
			wbOutput.write(os);
		}
	}
	
	/**
	 * Create a row in spreadsheet
	 * @param row Row number
	 * @param objects Objects to write in excel file
	 */
	private static void writeCells(Row row, Iterator<Cell> cellIterator) {
		int i = 0;
		while(cellIterator.hasNext()) {
			Cell prev = cellIterator.next();
			Cell cell = row.createCell(i, prev.getCellType());
			switch (prev.getCellType()) {
            case NUMERIC:
                cell.setCellValue(prev.getNumericCellValue());
                break;
            case STRING:
                cell.setCellValue(prev.getStringCellValue());
                break;
			default:
				cell.setCellValue("error");
				break;
			}
			i++;
		}
	}
}
