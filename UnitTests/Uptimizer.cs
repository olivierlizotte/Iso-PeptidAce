using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PeptidAce.Iso.Methods;

namespace PeptidAce.Iso.UnitTests
{
    public class PisTest
    {
        public static void Uptimize()
        {
            string fastaFile = @"C:\_IRIC\Data\NB\peptide.fasta";
            DBOptions options = PositionnalIsomerSolver.CreateOptions(fastaFile, @"C:\_IRIC\Data\NB\Units\", 8, 0.05, new PeptidAce.Utilities.Interfaces.ConSolCommandLine());
            
            options.dProduct =0.0886869235377232;
            options.dPrecursor =0.714634842572098;
            options.dMatchingProductFraction =0.432872176371921;
            options.dMatchingProduct =0.00492531899592156;
            options.dIntensityFraction =0.73908941342453;
            options.dIntensity = 0;// 0.687398171372431;
            options.dProtein =0.574124578188231;
            options.dPeptideScore =0.315866923572434;
            options.dFragmentScore = 0.0322849750669137;//*/



            string project = @"C:\_IRIC\Data\NB\ProjectTest_AllAce_Spiked_19Oct.csv";
            Samples samples = new Samples(project, 0, options);
            Uptimizer.Run(samples, options);
        }
    }

    public class Uptimizer
    {
        public static void Run(Samples samples, DBOptions options)
        {
            Ace ace = new Ace(options, samples);
            ace.Preload(false, false);            
            ace.PrepareQueries();

            bool keepUptimizing = true;
            int bestCount = 0;
            while (keepUptimizing)
            {
                try
                {
                    foreach (Query query in ace.AllQueries)
                    {
                        query.psms.Clear();
                        query.precursor.psms_AllPossibilities.Clear();
                        if (query.precursor.psms != null)
                            query.precursor.psms.Clear();
                    }//*/

                    Result rez = ace.LaunchSearch(ace.AllQueries);
                    int nbCorrectMatches = 0;

                    int nbTarget = 0;
                    foreach (Query query in rez.queries)
                        if (query.psms.Count > 0)
                            if (query.Target)
                            {
                                nbTarget++;
                                if (query.psms[0].Peptide.Sequence.CompareTo(query.sample.nameColumn) == 0)
                                    nbCorrectMatches++;
                            }
                    options.ConSole.WriteLine("NbCorrectMatches : " + nbCorrectMatches + "    (" + nbTarget + " targets)");//332
                    if (nbCorrectMatches > bestCount)
                    {
                        bestCount = nbCorrectMatches;
                        options.Save(options.OutputFolder + "Options_" + bestCount + ".csv");
                    }
                    options.RandomizeParams();
                }
                catch(Exception ex)
                {
                    Console.WriteLine(ex.Message);
                }
            }

            options.Save(options.OutputFolder + "UptimizedOptions.csv");
        }
    }
}
