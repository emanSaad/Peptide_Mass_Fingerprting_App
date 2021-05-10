/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfambiguity;
import java.util.ArrayList;
import org.biojava.bio.symbol.Symbol;


/**
 *
 * @author Eman
 */
public abstract class AmbiguousSymbol {

    String clearName; // the human readable name
    char charcode; // the character
    String bioName; // the biology name
    double averageMass; // the predefined average mass when the occurance of the symbol > 1
    int positionInSequence; // the position of the symbol in the sequence
    int lastUsedCorrespondenceChar; // this is used to keep track of during the replacement task
    ArrayList<Symbol> correspondingChars; // the set of symbols to replace this ambigouse one

   
    /**
     * get the name of symbol as biologists know it
     * @return
     */
    public String getClearName() {
        return clearName;
    }

    /**
     * To check if there still corresponding symbols for ambiguous symbol
     * @return
     */
    public boolean closed() {
        if (correspondingChars == null) {
            System.err.println("Ambigiuse symbol " + this + " Has no correspondece ");
            return true;
        }
        return (lastUsedCorrespondenceChar >= correspondingChars.size()-1);
    }
    

    /**
     * get the current corresponding symbol for the ambiguous symbol
     * @return
     */
    public Symbol currentSymbol() {
        
        return correspondingChars.get(lastUsedCorrespondenceChar);

    }

    /**
     * check next position of corresponding symbol 
     * if the sequence contains more than one ambiguous symbol
     * @return
     */
    public boolean nextPose() {
        if (closed()) {
            return false; // can not move 
        } else {
            lastUsedCorrespondenceChar++;
        }
        return true;
    }

    /**
     *Reset after finish replacing 
     */
    public void reset() {
        lastUsedCorrespondenceChar = 0;
    }

    /**
     * set clear name for ambiguous symbols
     * @param clearName
     */
    public void setClearName(String clearName) {
        this.clearName = clearName;
    }

    /**
     * Get the one char abbreviation for the ambiguous symbol 
     * @return
     */
    public char getCharcode() {
        return charcode;
    }

    /**
     * Set the one char abbreviation for the ambiguous symbol 
     * @param charcode
     */
    public void setCharcode(char charcode) {
        this.charcode = charcode;
    }

    /**
     * Get the biological name of ambiguous symbol, which is three letters abbreviation
     * @return
     */
    public String getBioName() {
        return bioName;
    }

    /**
     * Set the biological name of ambiguous symbol, which is three letters abbreviation
     * @param bioName
     */
    public void setBioName(String bioName) {
        this.bioName = bioName;
    }

    /**
     * Get the average mass of corresponding symbols of an ambiguous symbol
     * @return
     */
    public double getAverageMass() {
        return averageMass;
    }

    /**
     * Set the average mass of corresponding symbols of an ambiguous symbol
     * @param averageMass
     */
    public void setAverageMass(double averageMass) {
        this.averageMass = averageMass;
    }

    /**
     * get the position of ambiguous symbols in a sequence
     * @return position
     */
    public int getPositionInSequence() {
        return positionInSequence;
    }

    /**
     * Set the position of ambiguous symbols in a sequence
     * @param positionInSequence
     */
    public void setPositionInSequence(int positionInSequence) {
        this.positionInSequence = positionInSequence;
    }

    
    /**
     * Get the corresponding symbols of an ambiguous symbol
     * @return list of the corresponding symbols
     */
    public ArrayList<Symbol> getCorrespondingChars() {
        return correspondingChars;
    }

    /**
     * Set the corresponding symbols of an ambiguous symbol
     * @param correspondingChars
     */
    public void setCorrespondingChars(ArrayList<Symbol> correspondingChars) {
        this.correspondingChars = correspondingChars;
    }

    @Override
    public String toString() {
        return "AmbiguousSymbol{" + "clearName=" + clearName + ", charcode=" + charcode + ", bioName=" + bioName + ", averageMass=" + averageMass + ", positionInSequence=" + positionInSequence + ", lastUsedCorrespondenceChar=" + lastUsedCorrespondenceChar + '}';
    }
}
