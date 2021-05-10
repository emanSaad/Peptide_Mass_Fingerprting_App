/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfcore;

import java.util.Iterator;
import org.biojava.bio.proteomics.MassCalc;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.symbol.SymbolPropertyTable;

/** Calculate the molecular weight of protein
 *
 * @author Eman
 */
public class MassCalculator extends MassCalc{

    public MassCalculator(String isotopicType, boolean MH_PLUS) {
        super(isotopicType, MH_PLUS);
    }
  
    
    /**
     * Get the protein molecular weight
     * @param proteinSeq
     * @return
     * @throws IllegalSymbolException
     */
    public   double getMolecularWeight2(SymbolList proteinSeq)
     throws IllegalSymbolException {
        double pepMass = 0.0;
        Symbol gap = AlphabetManager.alphabetForName("PROTEIN").getGapSymbol();
        SymbolPropertyTable sPT =
            ProteinTools.getSymbolPropertyTable(SymbolPropertyTable.AVG_MASS);

        for (Iterator it = proteinSeq.iterator(); it.hasNext(); ) {
            Symbol s = (Symbol) it.next();
            if( s instanceof AtomicSymbol) {
                if(s==ProteinTools.o()) {
                    pepMass+=255.31;
                }
                else
                    if(s==ProteinTools.u()) {
                    pepMass+=168.05;
                }
                else {
                    pepMass += sPT.getDoubleValue(s);
                }
            } else {
                FiniteAlphabet matches = (FiniteAlphabet) s.getMatches();
                if(matches.size() == 0) {
                    if(s.equals(gap)) {
                        continue;
                    }
                    throw new IllegalSymbolException(
                        "Symbol " + s.getName() + " has no mass associated");
                } else {
                    int count = 0;
                    double mass = 0.0;
                    for(Iterator i = matches.iterator(); i.hasNext(); ) {
                        AtomicSymbol as = (AtomicSymbol) i.next();
                        //SELENOCYSTEINE and PYRROLYSINE
                        if(as==ProteinTools.o()) {
                            pepMass+=255.31;
                        }
                        else
                            if(as==ProteinTools.u()) {
                            pepMass+=168.05;
                        }
                        else {
                            mass += sPT.getDoubleValue(as);
                        }
                        count++;
                    }
                    pepMass += mass/((double)count);
                }
            }
        }

        //Calculate hydroxyl mass
        double termMass = calcTermMass(SymbolPropertyTable.AVG_MASS, false);

        if (pepMass != 0.0) {
            pepMass += termMass;
        }

        return pepMass;
     }
  
  private static double calcTermMass(String isotopicType, boolean MH_PLUS) {
        double termMass = 0.0;
        if (isotopicType.equals(SymbolPropertyTable.AVG_MASS)) {
            //Add the C-terminal OH and N-Term H
            termMass += Havg + Oavg + Havg;
            //Add the extra H
            if (MH_PLUS) {
                termMass += Havg;
            }
        }
        else if (isotopicType.equals(SymbolPropertyTable.MONO_MASS)) {
            //Add the C-terminal OH and N-Term H
            termMass += Hmono + Omono + Hmono;
            //Add the extra H
            if (MH_PLUS) {
                termMass += Hmono;
            }
        }
        return termMass;
    }  
}
