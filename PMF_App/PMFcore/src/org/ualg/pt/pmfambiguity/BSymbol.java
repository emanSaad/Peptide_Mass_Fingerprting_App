/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ualg.pt.pmfambiguity;

import java.util.ArrayList;
import org.ualg.pt.pmfcore.DigesterAndMassCalculater;

/**
 * Information of ambiguous symbol B
 * @author Eman
 */
public class BSymbol extends AmbiguousSymbol {

    public BSymbol(int index) {
        super();
        clearName = " symbol";
        charcode = 'B';
        averageMass = 132.6108;
        correspondingChars = new ArrayList(DigesterAndMassCalculater.getBsymbols());
        this.positionInSequence=index;
    }
}
