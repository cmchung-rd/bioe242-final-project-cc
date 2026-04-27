#ifndef TESTEPITHELIALFILLSIMULATION_HPP_
#define TESTEPITHELIALFILLSIMULATION_HPP_

/*
 * Direct C++ port of experiment7.6.py
 *
 * Epithelial tissue growth in a 60x8 rectangular domain using:
 *   - Vertex-based cell population (honeycomb mesh)
 *   - ContactInhibitionCellCycleModel
 *   - Neighbor-based proliferation inhibition (custom modifier)
 *   - NagaiHondaForce
 *   - Plane boundary conditions (rectangular box)
 *   - Incremental solve with coverage-based early stopping
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
#include "NagaiHondaForce.hpp"

// Boundary conditions
#include "PlaneBoundaryCondition.hpp"

// Modifiers
#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"

// Custom modifier
#include "NeighborDependentProliferationModifier.hpp"

// Utilities
#include "RandomNumberGenerator.hpp"

// Must be included last
#include "FakePetscSetup.hpp"

class TestEpithelialFillSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestEpithelialGrowthInRectangle()
    {
        // =========================================================
        //  Biological Parameters
        // =========================================================
        double g1Duration = 10.0;   // hours
        double sDuration  = 4.0;
        double g2Duration = 8.0;
        double mDuration  = 2.0;
        double totalCycleTime = g1Duration + sDuration + g2Duration + mDuration;

        // Contact inhibition
        double equilibriumVolume       = 1.0;
        double quiescentVolumeFraction = 0.60;

        // Neighbor-based inhibition
        unsigned neighborThreshold = 6;

        // =========================================================
        //  RNG seed (fixed for reproducibility)
        // =========================================================
        RandomNumberGenerator::Instance()->Reseed(123456789);

        // =========================================================
        //  Domain
        // =========================================================
        double boxW = 60.0;
        double boxH = 8.0;

        // =========================================================
        //  1. Generate mesh  (3 cols x 2 rows = 6 initial cells)
        // =========================================================
        HoneycombVertexMeshGenerator generator(3, 2);
        boost::shared_ptr<MutableVertexMesh<2, 2> > p_mesh = generator.GetMesh();

        // Mesh rearrangement thresholds for stability at large cell counts
        p_mesh->SetCellRearrangementThreshold(0.05);
        p_mesh->SetT2Threshold(0.005);
        p_mesh->SetCheckForT3Swaps(false);

        // ----------------------------------------------------------
        //  Position initial cells: left-aligned, vertically centred
        // ----------------------------------------------------------
        c_vector<double, 2> mesh_min;
        c_vector<double, 2> mesh_max;
        mesh_min[0] =  1e30; mesh_min[1] =  1e30;
        mesh_max[0] = -1e30; mesh_max[1] = -1e30;

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

        // Target: left edge of mesh at x=0, centred at y = boxH/2
        c_vector<double, 2> target;
        target[0] = mesh_center[0] - mesh_min[0];   // shift so leftmost node → x=0
        target[1] = boxH / 2.0;

        c_vector<double, 2> offset = target - mesh_center;

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2> new_loc = p_mesh->GetNode(i)->rGetLocation() + offset;
            ChastePoint<2> point(new_loc);
            p_mesh->GetNode(i)->SetPoint(point);
        }

        // =========================================================
        //  2. Create cells with ContactInhibitionCellCycleModel
        // =========================================================
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        CellsGenerator<ContactInhibitionCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        // Configure each cell's cycle model
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

            // Stagger birth times so divisions aren't synchronised
            double birthTime = -totalCycleTime * RandomNumberGenerator::Instance()->ranf();
            cells[i]->SetBirthTime(birthTime);
        }

        // =========================================================
        //  3. Cell population
        // =========================================================
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // =========================================================
        //  4. Simulation
        // =========================================================
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("EpithelialRheology/FillSimulation");
        simulator.SetSamplingTimestepMultiple(2000);   // output every 2000×dt = 1 h
        simulator.SetDt(0.0005);

        // -- Modifiers (order matters: volume tracking first) ------
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetReferenceTargetArea(equilibriumVolume);
        simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR(NeighborDependentProliferationModifier<2>, p_neighbor_modifier);
        p_neighbor_modifier->SetNeighborThreshold(neighborThreshold);
        simulator.AddSimulationModifier(p_neighbor_modifier);

        // -- Force --------------------------------------------------
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(0.5);
        simulator.AddForce(p_force);

        // -- Boundary conditions (rectangular box) ------------------
        typedef PlaneBoundaryCondition<2, 2> PlaneBC;

        // Bottom wall  y = 0
        c_vector<double, 2> lower_point = zero_vector<double>(2);
        c_vector<double, 2> normal_down = zero_vector<double>(2);
        normal_down[1] = -1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_lower,
                       (&cell_population, lower_point, normal_down));
        simulator.AddCellPopulationBoundaryCondition(p_bc_lower);

        // Top wall  y = boxH
        c_vector<double, 2> upper_point = zero_vector<double>(2);
        upper_point[1] = boxH;
        c_vector<double, 2> normal_up = zero_vector<double>(2);
        normal_up[1] = 1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_upper,
                       (&cell_population, upper_point, normal_up));
        simulator.AddCellPopulationBoundaryCondition(p_bc_upper);

        // Left wall  x = 0
        c_vector<double, 2> left_point = zero_vector<double>(2);
        c_vector<double, 2> normal_left = zero_vector<double>(2);
        normal_left[0] = -1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_left,
                       (&cell_population, left_point, normal_left));
        simulator.AddCellPopulationBoundaryCondition(p_bc_left);

        // Right wall  x = boxW
        c_vector<double, 2> right_point = zero_vector<double>(2);
        right_point[0] = boxW;
        c_vector<double, 2> normal_right = zero_vector<double>(2);
        normal_right[0] = 1.0;
        MAKE_PTR_ARGS(PlaneBC, p_bc_right,
                       (&cell_population, right_point, normal_right));
        simulator.AddCellPopulationBoundaryCondition(p_bc_right);

        // =========================================================
        //  5. Incremental solve with coverage-based early stopping
        // =========================================================
        double targetBoxArea = boxW * boxH;
        double timeIncrement = 1.0;    // check every 1 hour
        double currentTime   = 0.0;
        double maxTime       = 600.0;

        std::cout << "Starting simulation (contact inhibition ON, no apoptosis)...\n";
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

            // Compute coverage
            double totalCellArea = 0.0;
            for (unsigned elem_idx = 0; elem_idx < p_mesh->GetNumElements(); elem_idx++)
            {
                totalCellArea += p_mesh->GetVolumeOfElement(elem_idx);
            }
            double coverage = totalCellArea / targetBoxArea;
            unsigned numCells = cell_population.GetNumRealCells();

            std::cout << "t=" << currentTime << " h | cells=" << numCells
                      << " | area=" << totalCellArea << "/" << targetBoxArea
                      << " | coverage=" << coverage * 100.0 << "%\n";

            if (coverage >= 1.0)
            {
                std::cout << "Box fully covered at t=" << currentTime
                          << " h. Stopping.\n";
                break;
            }
        }

        if (currentTime >= maxTime)
        {
            std::cout << "Reached max time (" << maxTime
                      << " h) without full coverage.\n";
        }

        std::cout << "Final cell count: " << cell_population.GetNumRealCells() << "\n";
    }
};

#endif // TESTEPITHELIALFILLSIMULATION_HPP_
