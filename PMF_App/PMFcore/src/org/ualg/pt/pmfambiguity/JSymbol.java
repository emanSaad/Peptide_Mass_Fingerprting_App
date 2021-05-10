/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfambiguity;

import java.util.ArrayList;
import org.ualg.pt.pmfcore.DigesterAndMassCalculater;

/**
 * Information of ambiguous symbol J
 * @see AmbiguousSymbol
 * @author Eman
 * 
 */
public class JSymbol extends AmbiguousSymbol {

    public JSymbol(int index) {
        super();
        clearName = "J symbol";
        charcode = 'J';
        averageMass = 131.1736;
        correspondingChars = new ArrayList(DigesterAndMassCalculater.getJsymbols());
        this.positionInSequence=index;
    }
}
