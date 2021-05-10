/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfambiguity;

import java.util.ArrayList;
import org.ualg.pt.pmfcore.DigesterAndMassCalculater;

/**
 * Information of ambiguous symbol X
 * @author Eman
 */
public class XSymbol extends AmbiguousSymbol {

    public XSymbol(int index) {
        super();
        clearName = "X symbol";
        charcode = 'X';
        averageMass = 110;
        correspondingChars = new ArrayList(DigesterAndMassCalculater.getXsymbols());
        this.positionInSequence=index;
    }
}
