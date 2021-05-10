/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfambiguity;

import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.utils.ChangeVetoException;
import org.ualg.pt.pmfcore.Protein;

/**
 *
 * @author Eman
 */
public class AmbiguityResolver {

    AmbiguousSymbol[] ambiguityMap;
    int nextIndex=0;
    /**
     * Create two arrays one for the symbols to be replaced by a set of symbols.
     * The other map is for those symbols which will be replaced by the average
     * mass
     *
     * @param protein the sequence to calculate it's ambiguity maps
     * @param unknown the map of unknown symbols and their occurrence in the
     * symbol
     */
    public void calcMaps(Protein protein, Map<Symbol, LinkedList<Integer>> unknown) {


        int singleSymbolCount = 0;
        int multipleSymbolCount = 0;

        // calcualte the number of multiple symbols
        for (Iterator<Map.Entry<Symbol, LinkedList<Integer>>> it = unknown.entrySet().iterator(); it.hasNext();) {

            Map.Entry<Symbol, LinkedList<Integer>> entry = it.next();

            if (entry.getValue().size() > 1) {
                multipleSymbolCount++;
            }
        }

        // create two arrays
        int length=unknown.size() - multipleSymbolCount;
        if(length>0) {
            ambiguityMap = new AmbiguousSymbol[unknown.size() - multipleSymbolCount];
        }
   
        // calculate the molar replacement map, i.e the average mass 
        // calculate the symbol replacement map

        int i = 0;
        for (Map.Entry<Symbol, LinkedList<Integer>> entry : unknown.entrySet()) {
            LinkedList<Integer> l = entry.getValue();           

            if (l.size() == 1) { // only one symbol 
                ambiguityMap[i++] = fromSymbol2AmbiguouseSymbol(entry.getKey(), entry.getValue().get(0)); // we are expecting one occurence.
            }
        }



    }

   
    /**
     * Get the ambiguous symbols
     * @return ambiguous symbols
     */
    public AmbiguousSymbol[] getAmbiguityMap() {
        return ambiguityMap;
    }

    
    /**
     * return the ambiguous symbols by their names
     * @return
     */
    public String getMapAsString(){
        String map="";
       for(int i=0;i<ambiguityMap.length;i++){
           map+=ambiguityMap[i].currentSymbol().getName();
       }  
       return map;
    }
    
    /**
     * Set the ambiguous symbols
     * @param ambiguityMap
     */
    public void setAmbiguityMap(AmbiguousSymbol[] ambiguityMap) {
        this.ambiguityMap = ambiguityMap;
    }
    

    /**
     * Check which symbol is ambiguous symbol, each one has more than one corresponding standard symbols
     * @return ambiguous symbol
     */
    private AmbiguousSymbol fromSymbol2AmbiguouseSymbol(Symbol s, int index) {
        try {
            
            XSymbol xSymbol = new XSymbol(index);
            if (s.getMatches().contains(ProteinTools.getAlphabet().getAmbiguity(new HashSet(xSymbol.getCorrespondingChars())))) {
                return xSymbol;
            }
            

            BSymbol bSymbol = new BSymbol(index);
            if (s.getMatches().contains(ProteinTools.getAlphabet().getAmbiguity(new HashSet(bSymbol.getCorrespondingChars())))) {
                return bSymbol;
            }

            JSymbol jSymbol = new JSymbol(index);
            if (s.getMatches().contains(ProteinTools.getAlphabet().getAmbiguity(new HashSet(jSymbol.getCorrespondingChars())))) {
                return jSymbol;
            }


            ZSymbol zSymbol = new ZSymbol(index);
            if (s.getMatches().contains(ProteinTools.getAlphabet().getAmbiguity(new HashSet(zSymbol.getCorrespondingChars())))) {
                return zSymbol;
            }

        } catch (IllegalSymbolException ex) {
            Logger.getLogger(AmbiguityResolver.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    
    
    /**
     * Reset the ambiguity map
     */
    public void resetMap() {
        for (int i = 0; i < ambiguityMap.length; i++) {
            ambiguityMap[i].reset();
        }

    }

    /**
     * resolve the sequence to generate a new one
     *
     * @param protein
     * @return true if success, false otherwise
     */
    public Protein resolveSequence(Protein protein) {

        Protein newProtein=null;
        try {
            newProtein = new Protein(protein);
            // if the last ambiguous symbol has been replaced return false;
        } catch (BioException ex) {
            Logger.getLogger(AmbiguityResolver.class.getName()).log(Level.SEVERE, null, ex);
        }

        if (ambiguityMap == null) {
            return null; // there is nothing to replace
        }          
        
        // construct the next replacement 
           for (int i = 0; i < ambiguityMap.length; i++) {
            try {
               
                Edit edRest = new Edit(ambiguityMap[i].getPositionInSequence(), protein.getSequence().getAlphabet(), ambiguityMap[i].currentSymbol());
                try {
                    newProtein.getSequence().edit(edRest);
                } catch (IndexOutOfBoundsException ex) {
                    Logger.getLogger(AmbiguityResolver.class.getName()).log(Level.SEVERE, null, ex);
                } catch (IllegalAlphabetException ex) {
                    Logger.getLogger(AmbiguityResolver.class.getName()).log(Level.SEVERE, null, ex);
                } catch (ChangeVetoException ex) {
                    Logger.getLogger(AmbiguityResolver.class.getName()).log(Level.SEVERE, null, ex);
                }
            } catch (IllegalSymbolException ex) {
                Logger.getLogger(AmbiguityResolver.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
           
// use this commented code to check if everything is working fine           
           
//         for (int i = 0; i < ambiguityMap.length ; i++)
//             System.out.print(ambiguityMap[i].lastUsedCorrespondenceChar +"\t" );
//
//              System.out.println();
                                      
          boolean allclosed=true;
        for (int i = 0; i < ambiguityMap.length; i++) {
         if (!ambiguityMap[i].closed()){
             allclosed=false;
             break;
         }
        }
        if( allclosed) {
            return null;
        } // no avilable symbols to replace.
                 
         for(int i=0;i<ambiguityMap.length;i++){
             if(ambiguityMap[i].closed()){
                 ambiguityMap[i].reset();
                 continue;                 
             }
         else
             {
                 ambiguityMap[i].nextPose();
                               
             }
            break;
         }   
        // generate a new sequence    
        return newProtein;
    }
}
