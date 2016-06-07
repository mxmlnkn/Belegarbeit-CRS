#include "ulf.h"
extern "C" {
    void postProcess_PV( equation *eq, field * rhs __attribute__((unused)) )
    {
        ASSERT_OUT("You are using post-processing function to calculate PV based on selected species");
        printf("You are using post-processing function to calculate PV based on selected species\n");
        parser * speciesPrsr    = EQ_PARSER ("species"); 
		ASSERT_MSG(speciesPrsr, "postProcessing_PV: 'species' required in 'postProcessing'!");
        problem & ulfMain = EQ_PRB;
        database * ulfDB = ulfMain.getDatabase();
        field *PV = ulfDB->getField("PV");
        ulfDB->setOptions("PV", WRITE_FIELD);
        mixture & fld  = EQ_MIX; 
        field **Yfld = fld.getY();
        int nSpecies = EQ_NSPEC;
        double alpha[nSpecies];
        for(int i=0;i<nSpecies;i++)
        {
            alpha[i] = 0.0;
        }
        forAllSeq(*speciesPrsr, i)
        {
            char *name   = (*speciesPrsr)[i].getKey();
            double value = (*speciesPrsr)[i].getDouble();
            int iSpec    = fld.speciesIndex(name);
            if(iSpec > -1)
            {
                alpha[iSpec] = value;
            }
        }
        forAllSeq(*PV,i)
        {
            (*PV)[i] = 0.0;
            for(int j=0;j<nSpecies;j++)
            {
                (*PV)[i] += alpha[j] * (*Yfld[j])[i];
            }
        }

    }
    // SDR ... skalare dissipationsrate
    void postProcess_SDRConfigure(equation *eq)
    {
        EQ_REG_F("Z");      // 0
        EQ_REG_F("F");      // 1 seems to be the Y->f<ields
        EQ_REG_F("rho");    // 2
        EQ_REG_F("lambda"); // 3
        EQ_REG_F("cpMean"); // 4
        EQ_REG_F("chi");    // 5
    }
    void postProcess_SDR(equation *eq, field *rhs __attribute__((unused)))
    {
        // get fields
        field &Z = *EQ_F(0);
        field &F = *EQ_F(1);
        field &rho = *EQ_F(2);
        field &lambda = *EQ_F(3);
        field &cpMean = *EQ_F(4);
        field &chi = *EQ_F(5);

        field dZdx = Z.derivative(0, FIRST, &F);
        chi = 2.0 * lambda / (rho * cpMean) * dZdx * dZdx; // Le_Z=1
    }
}

