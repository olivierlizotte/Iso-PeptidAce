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
            DBOptions dbOptions = PositionnalIsomerSolver.CreateOptions(fastaFile, @"C:\_IRIC\Data\NB\Units2\", 8, 0.05, new PeptidAce.Utilities.Interfaces.ConSolCommandLine());

            dbOptions.dProduct = 0.0917981081138356;
            dbOptions.dPrecursor = 0.345789190542786;
            dbOptions.dMatchingProductFraction = 0.427418045898628;
            dbOptions.dMatchingProduct = 0;
            dbOptions.dIntensityFraction = 0.429418127252449;
            dbOptions.dIntensity = 0;
            dbOptions.dProtein = 0.692270441303156;
            dbOptions.dPeptideScore = 0.636739763262095;
            dbOptions.dFragmentScore = 0.0229058195943506;
            
            dbOptions.fullFragment = new FullFragments(false, false, false);
            string project = @"C:\_IRIC\Data\NB\ProjectTest_AllAce_Spiked_QEPlus_Apr21.csv";
            Samples samples = new Samples(project, 0, dbOptions);
            Uptimizer.Run(samples, dbOptions);

            dbOptions.Save(dbOptions.OutputFolder + "UptimizedOptions.csv");
        }
    }

    public class Uptimizer
    {
        public static void Run(Samples samples, DBOptions options)
        {
            //PeptidAce.Utilities.Methods.UptimizeOptions upper = new Utilities.Methods.UptimizeOptions(options, samples);
            //upper.Run();
            
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
                    }

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
            //*/
        }
    }
}
