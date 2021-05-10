/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import org.biojava.bio.BioException;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojavax.SimpleNamespace;
import org.biojavax.bio.seq.RichSequence;
import org.biojavax.bio.seq.RichSequenceIterator;


/**
 *
 * @author Eman
 */
public  abstract class DataReader {
    
    //File massesFile;
    // Read the masses file and return arraylist masses 
    /**
     * Read experimental masses file that contains double values 
     * @param massesFile
     * @return Double list of masses 
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static ArrayList<Double> readMasses(File massesFile) throws FileNotFoundException, IOException
    {
        ArrayList <Double> massesArrayList;       
        BufferedReader reader = new BufferedReader(new FileReader(massesFile));
        String line;
        massesArrayList= new ArrayList <Double>();
        while((line = reader.readLine()) != null) 
        {
            if ("".equals(line)) {
                 System.out.println("You select invalid masses file, please check if it contains empty lines.....");
                 System.exit(0);
            }
            else
            {
                massesArrayList.add(Double.valueOf(line));     
            }              
        }
        System.out.println("Masses File: "+ massesFile.getAbsolutePath());
       //System.out.println("Eperimental Data: \t"+massesArrayList);
       return massesArrayList; 
    }
    
    
    /**
     
     * Read protein sequences from FASTA files and return a sequence
     * @param file
     * @return a list of sequences found in the file
     */
    public static ArrayList <RichSequence> readSequence(File file)
    {
        RichSequence rec = null;
        BufferedReader br =null;
        ArrayList<RichSequence> sequences= new ArrayList<RichSequence>();
        
        try {                       
            br = new BufferedReader(new FileReader(file));
            Alphabet alpha = AlphabetManager.alphabetForName("PROTEIN");
            SimpleNamespace ns = new SimpleNamespace("biojava");
            
            RichSequenceIterator iterator = RichSequence.IOTools.readFasta(br,alpha.getTokenization("token"), ns);
            while (iterator.hasNext()) 
            {
                rec = iterator.nextRichSequence();
                sequences.add(rec);
                //System.out.println(rec.getName());
                //System.out.println(rec.length());
            }
          
        } catch (FileNotFoundException ex) {
           // Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("File not found!");
        } catch (BioException ex) {
            //Logger.getLogger(DataReader.class.getName()).log(Level.SEVERE, null, ex);
            System.err.println("Can not load sequence DataReader Class:"+ ex.getMessage());
        }
         return sequences;
    }    
}
    
