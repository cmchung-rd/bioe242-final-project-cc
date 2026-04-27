#ifndef TESTPROLIFERATION_HPP_
#define TESTPROLIFERATION_HPP_

/*
 * Proliferation parameter sweep — 30x4 domain
 *
 * Identical to TestControl_30x4 except total cell cycle time is scaled uniformly.
 * All other parameters are held at control values.
 *
 * Control cycle: g1=10, S=4, G2=8, M=2 (total=24 h)
 *
 * Variants — scale factor applied to all phase durations:
 *   x0.10, x0.20, x0.33, x0.50, x0.75  <- faster proliferation
 *   x1.10, x1.25, x1.50, x1.75, x2.00  <- slower proliferation
 *
 * Output: EpithelialRheology/Proliferation/CycleScale_<factor>/
 */

#include <cxxtest/TestSuite.h>

// Chaste core
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"

// Mesh
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "ChastePoint.hpp"

// Cell population
#include "VertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"

// Cell cycle
#include "ContactInhibitionCellCycleModel.hpp"

// Cell types
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

// Simulation
#include "OffLatticeSimulation.hpp"

// Forces
#include "FarhadifarForce.hpp"
#include "AbstractForce.hpp"

// Boundary conditions
#include "PlaneBoundaryCondition.hpp"

// Cell killers
#include "RandomCellKiller.hpp"
#include "T2SwapCellKiller.hpp"

// Modifiers
#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"

// Custom modifier
#include "NeighborDependentProliferationModifier.hpp"

// Output writers
#include "CellVolumesWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"

// Utilities
#include "RandomNumberGenerator.hpp"

// Must be last
#include "FakePetscSetup.hpp"

// =========================================================================
// CUSTOM FORCE: Active Motility with Wall Drag (Parabolic Flow Profile)
// =========================================================================
template<unsigned DIM>
class ParabolicLeadingEdgeForce : public AbstractForce<DIM>
{
private:
    double mActiveForceMagnitude;
    double mBoxHeight;

public:
    ParabolicLeadingEdgeForce(double activeForceMagnitude, double boxHeight)
        : AbstractForce<DIM>(),
          mActiveForceMagnitude(activeForceMagnitude),
          mBoxHeight(boxHeight)
    {}

    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        double max_x = -1e30;
        for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            if (node_iter->rGetLocation()[0] > max_x)
                max_x = node_iter->rGetLocation()[0];
        }

        for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            double x = node_iter->rGetLocation()[0];
            double y = node_iter->rGetLocation()[1];

            if (x > max_x - 1.5)
            {
                double center_y = mBoxHeight / 2.0;
                double normalized_dist = (y - center_y) / center_y;
                double scale = 1.0 - (normalized_dist * normalized_dist);
                if (scale < 0.0) scale = 0.0;

                c_vector<double, DIM> force = zero_vector<double>(DIM);
                force[0] = mActiveForceMagnitude * scale;
                node_iter->AddAppliedForceContribution(force);
            }
        }
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<ActiveForceMagnitude>" << mActiveForceMagnitude << "</ActiveForceMagnitude>\n";
        *rParamsFile << "\t\t\t<BoxHeight>" << mBoxHeight << "</BoxHeight>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};


class TestProliferation : public AbstractCellBasedTestSuite
{
private:

    /*
     * Runs the 30x4 fill simulation with cycle durations scaled by cycleScaleFactor.
     * All other parameters are fixed at control values.
     *
     * @param cycleScaleFactor  multiplier applied to all cell cycle phase durations
     * @param outputLabel       subfolder name, e.g. "CycleScale_0.10"
     */
    void RunProliferationVariant(double cycleScaleFactor,
                                 const std::string& outputLabel,
                                 unsigned seed)
    {
        // =========================================================
        //  Parameters (all fixed at control values except cycle durations)
        // =========================================================

        // Domain
        const double boxW = 15.0;
        const double boxH =  4.0;

        // FarhadifarForce (control values)
        const double areaElasticity         = 10.0;
        const double perimeterContractility = 0.04;
        const double lineTension            = 0.12;
        const double boundaryLineTension    = 0.12;

        // Active Motility Parameter
        const double activeMotilityForceMagnitude = 0.1;

        // Cell cycle — base durations scaled by cycleScaleFactor
        const double g1Duration  = 10.0 * cycleScaleFactor;
        const double sDuration   =  4.0 * cycleScaleFactor;
        const double g2Duration  =  8.0 * cycleScaleFactor;
        const double mDuration   =  2.0 * cycleScaleFactor;
        const double totalCycleTime = g1Duration + sDuration + g2Duration + mDuration;
        const double equilibriumVolume       = 1.0;
        const double quiescentVolumeFraction = 0.10;

        // Neighbor-based proliferation arrest
        const unsigned neighborThreshold = 6;

        // Cell death (disabled)
        const double deathProbabilityInAnHour = 0.0;

        // Simulation
        const double   dt               = 0.0005;
        const unsigned samplingMultiple = 2000;
        const double   timeIncrement    = 1.0;
        const double   maxTime          = 1200.0;

        // =========================================================
        //  RNG seed (fixed for reproducibility)
        // =========================================================
        RandomNumberGenerator::Instance()->Reseed(seed);

        // =========================================================
        //  1. Generate mesh  (3 cols x 2 rows = 6 initial cells)
        // =========================================================
        HoneycombVertexMeshGenerator generator(3, 2);
        boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = generator.GetMesh();

        p_mesh->SetCellRearrangementThreshold(0.01);
        p_mesh->SetT2Threshold(0.001);
        p_mesh->SetCheckForT3Swaps(false);

        // Position initial cells: left-aligned, vertically centred in the box
        c_vector<double, 2> mesh_min, mesh_max;
        mesh_min[0] =  1e30;  mesh_min[1] =  1e30;
        mesh_max[0] = -1e30;  mesh_max[1] = -1e30;

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2> loc = p_mesh->GetNode(i)->rGetLocation();
            for (unsigned d = 0; d < 2; d++)
            {
                if (loc[d] < mesh_min[d]) mesh_min[d] = loc[d];
                if (loc[d] > mesh_max[d]) mesh_max[d] = loc[d];
            }
        }

        c_vector<double, 2> mesh_center;
        mesh_center[0] = 0.5 * (mesh_min[0] + mesh_max[0]);
        mesh_center[1] = 0.5 * (mesh_min[1] + mesh_max[1]);

        c_vector<double, 2> target;
        target[0] = mesh_center[0] - mesh_min[0];
        target[1] = boxH / 2.0;

        c_vector<double, 2> offset = target - mesh_center;

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2> new_loc = p_mesh->GetNode(i)->rGetLocation() + offset;
            ChastePoint<2> point(new_loc);
            p_mesh->GetNode(i)->SetPoint(point);
        }

        // =========================================================
        //  2. Create cells
        // =========================================================
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        CellsGenerator<ContactInhibitionCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        for (unsigned i = 0; i < cells.size(); i++)
        {
            ContactInhibitionCellCycleModel* p_model =
                static_cast<ContactInhibitionCellCycleModel*>(cells[i]->GetCellCycleModel());

            p_model->SetTransitCellG1Duration(g1Duration);
            p_model->SetSDuration(sDuration);
            p_model->SetG2Duration(g2Duration);
            p_model->SetMDuration(mDuration);
            p_model->SetEquilibriumVolume(equilibriumVolume);
            p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);

            double birthTime = -totalCycleTime * RandomNumberGenerator::Instance()->ranf();
            cells[i]->SetBirthTime(birthTime);
        }

        // =========================================================
        //  3. Cell population
        // =========================================================
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddPopulationWriter<NodeVelocityWriter>();
        cell_population.AddPopulationWriter<VertexT1SwapLocationsWriter>();

        // =========================================================
        //  4. Simulation setup
        // =========================================================
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EpithelialRheology/Proliferation/" + outputLabel);
        simulator.SetSamplingTimestepMultiple(samplingMultiple);
        simulator.SetDt(dt);

        // -- Modifiers --
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetReferenceTargetArea(equilibriumVolume);
        simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR(NeighborDependentProliferationModifier<2>, p_neighbor_modifier);
        p_neighbor_modifier->SetNeighborThreshold(neighborThreshold);
        simulator.AddSimulationModifier(p_neighbor_modifier);

        // -- Force --
        MAKE_PTR(FarhadifarForce<2>, p_force);
        p_force->SetAreaElasticityParameter(areaElasticity);
        p_force->SetPerimeterContractilityParameter(perimeterContractility);
        p_force->SetLineTensionParameter(lineTension);
        p_force->SetBoundaryLineTensionParameter(boundaryLineTension);
        simulator.AddForce(p_force);

        MAKE_PTR_ARGS(ParabolicLeadingEdgeForce<2>, p_motility_force, (activeMotilityForceMagnitude, boxH));
        simulator.AddForce(p_motility_force);

        // -- Boundary conditions --
        typedef PlaneBoundaryCondition<2, 2> PlaneBC;

        c_vector<double, 2> lower_point = zero_vector<double>(2);
        c_vector<double, 2> normal_down = zero_vector<double>(2);
        normal_down[1] = -1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_lower, (&cell_population, lower_point, normal_down));
        simulator.AddCellPopulationBoundaryCondition(p_bc_lower);

        c_vector<double, 2> upper_point = zero_vector<double>(2);
        upper_point[1] = boxH;
        c_vector<double, 2> normal_up = zero_vector<double>(2);
        normal_up[1] = 1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_upper, (&cell_population, upper_point, normal_up));
        simulator.AddCellPopulationBoundaryCondition(p_bc_upper);

        c_vector<double, 2> left_point = zero_vector<double>(2);
        c_vector<double, 2> normal_left = zero_vector<double>(2);
        normal_left[0] = -1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_left, (&cell_population, left_point, normal_left));
        simulator.AddCellPopulationBoundaryCondition(p_bc_left);

        c_vector<double, 2> right_point = zero_vector<double>(2);
        right_point[0] = boxW;
        c_vector<double, 2> normal_right = zero_vector<double>(2);
        normal_right[0] = 1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_right, (&cell_population, right_point, normal_right));
        simulator.AddCellPopulationBoundaryCondition(p_bc_right);

        // -- Cell killers --
        MAKE_PTR_ARGS(RandomCellKiller<2>, p_random_killer,
                      (&cell_population, deathProbabilityInAnHour));
        simulator.AddCellKiller(p_random_killer);

        MAKE_PTR_ARGS(T2SwapCellKiller<2>, p_t2_killer, (&cell_population));
        simulator.AddCellKiller(p_t2_killer);

        // =========================================================
        //  5. Incremental solve — stops at 100% coverage of 120 units^2
        // =========================================================
        const double targetBoxArea = boxW * boxH;   // 120.0 units^2
        double currentTime = 0.0;

        std::cout << "Starting simulation: " << outputLabel << "\n";
        std::cout << "cycleScaleFactor = " << cycleScaleFactor
                  << " (total cycle = " << totalCycleTime << " h)\n";
        std::cout << "Target area: " << targetBoxArea << " units^2\n";
        std::cout << "Initial cell count: " << cell_population.GetNumRealCells() << "\n";

        while (currentTime < maxTime)
        {
            double nextTime = std::min(currentTime + timeIncrement, maxTime);
            simulator.SetEndTime(nextTime);

            try
            {
                simulator.Solve();
            }
            catch (Exception& e)
            {
                std::cout << "\nMesh error at t=" << nextTime
                          << " h: " << e.GetMessage() << "\n";
                std::cout << "Stopping early due to mesh topology failure.\n";
                break;
            }

            currentTime = nextTime;

            double totalCellArea = 0.0;
            for (unsigned elem_idx = 0; elem_idx < p_mesh->GetNumElements(); elem_idx++)
            {
                totalCellArea += p_mesh->GetVolumeOfElement(elem_idx);
            }
            double coverage = totalCellArea / targetBoxArea;

            std::cout << "t=" << currentTime
                      << " h | cells=" << cell_population.GetNumRealCells()
                      << " | area=" << totalCellArea << "/" << targetBoxArea
                      << " | coverage=" << coverage * 100.0 << "%\n";

            if (coverage >= 1.0)
            {
                std::cout << "Box fully covered (100%) at t=" << currentTime
                          << " h. Stopping.\n";
                break;
            }
        }

        if (currentTime >= maxTime)
        {
            std::cout << "Reached safety cap (" << maxTime
                      << " h) without full coverage.\n";
        }

        std::cout << "Final cell count: " << cell_population.GetNumRealCells() << "\n";
    }

public:

    // ── 5 faster variants — 3 replicates each ─────────────────────────

    void TestCycleScale_0p10_rep1() { RunProliferationVariant(0.10, "CycleScale_0.10_rep1", 1001); }
    void TestCycleScale_0p10_rep2() { RunProliferationVariant(0.10, "CycleScale_0.10_rep2", 2002); }
    void TestCycleScale_0p10_rep3() { RunProliferationVariant(0.10, "CycleScale_0.10_rep3", 3003); }

    void TestCycleScale_0p20_rep1() { RunProliferationVariant(0.20, "CycleScale_0.20_rep1", 1001); }
    void TestCycleScale_0p20_rep2() { RunProliferationVariant(0.20, "CycleScale_0.20_rep2", 2002); }
    void TestCycleScale_0p20_rep3() { RunProliferationVariant(0.20, "CycleScale_0.20_rep3", 3003); }

    void TestCycleScale_0p33_rep1() { RunProliferationVariant(0.33, "CycleScale_0.33_rep1", 1001); }
    void TestCycleScale_0p33_rep2() { RunProliferationVariant(0.33, "CycleScale_0.33_rep2", 2002); }
    void TestCycleScale_0p33_rep3() { RunProliferationVariant(0.33, "CycleScale_0.33_rep3", 3003); }

    void TestCycleScale_0p50_rep1() { RunProliferationVariant(0.50, "CycleScale_0.50_rep1", 1001); }
    void TestCycleScale_0p50_rep2() { RunProliferationVariant(0.50, "CycleScale_0.50_rep2", 2002); }
    void TestCycleScale_0p50_rep3() { RunProliferationVariant(0.50, "CycleScale_0.50_rep3", 3003); }

    void TestCycleScale_0p75_rep1() { RunProliferationVariant(0.75, "CycleScale_0.75_rep1", 1001); }
    void TestCycleScale_0p75_rep2() { RunProliferationVariant(0.75, "CycleScale_0.75_rep2", 2002); }
    void TestCycleScale_0p75_rep3() { RunProliferationVariant(0.75, "CycleScale_0.75_rep3", 3003); }

    // ── 5 slower variants — 3 replicates each ─────────────────────────

    void TestCycleScale_1p10_rep1() { RunProliferationVariant(1.10, "CycleScale_1.10_rep1", 1001); }
    void TestCycleScale_1p10_rep2() { RunProliferationVariant(1.10, "CycleScale_1.10_rep2", 2002); }
    void TestCycleScale_1p10_rep3() { RunProliferationVariant(1.10, "CycleScale_1.10_rep3", 3003); }

    void TestCycleScale_1p25_rep1() { RunProliferationVariant(1.25, "CycleScale_1.25_rep1", 1001); }
    void TestCycleScale_1p25_rep2() { RunProliferationVariant(1.25, "CycleScale_1.25_rep2", 2002); }
    void TestCycleScale_1p25_rep3() { RunProliferationVariant(1.25, "CycleScale_1.25_rep3", 3003); }

    void TestCycleScale_1p50_rep1() { RunProliferationVariant(1.50, "CycleScale_1.50_rep1", 1001); }
    void TestCycleScale_1p50_rep2() { RunProliferationVariant(1.50, "CycleScale_1.50_rep2", 2002); }
    void TestCycleScale_1p50_rep3() { RunProliferationVariant(1.50, "CycleScale_1.50_rep3", 3003); }

    void TestCycleScale_1p75_rep1() { RunProliferationVariant(1.75, "CycleScale_1.75_rep1", 1001); }
    void TestCycleScale_1p75_rep2() { RunProliferationVariant(1.75, "CycleScale_1.75_rep2", 2002); }
    void TestCycleScale_1p75_rep3() { RunProliferationVariant(1.75, "CycleScale_1.75_rep3", 3003); }

    void TestCycleScale_2p00_rep1() { RunProliferationVariant(2.00, "CycleScale_2.00_rep1", 1001); }
    void TestCycleScale_2p00_rep2() { RunProliferationVariant(2.00, "CycleScale_2.00_rep2", 2002); }
    void TestCycleScale_2p00_rep3() { RunProliferationVariant(2.00, "CycleScale_2.00_rep3", 3003); }
};

#endif // TESTPROLIFERATION_HPP_
