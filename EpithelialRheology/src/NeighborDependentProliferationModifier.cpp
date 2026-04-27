#include "NeighborDependentProliferationModifier.hpp"
#include "SmartPointers.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MutableVertexMesh.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template<unsigned DIM>
NeighborDependentProliferationModifier<DIM>::NeighborDependentProliferationModifier()
    : AbstractCellBasedSimulationModifier<DIM, DIM>(),
      mNeighborThreshold(6)
{
}

template<unsigned DIM>
void NeighborDependentProliferationModifier<DIM>::SetNeighborThreshold(unsigned threshold)
{
    mNeighborThreshold = threshold;
}

template<unsigned DIM>
unsigned NeighborDependentProliferationModifier<DIM>::GetNeighborThreshold() const
{
    return mNeighborThreshold;
}

template<unsigned DIM>
void NeighborDependentProliferationModifier<DIM>::UpdateAtEndOfTimeStep(
    AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    // This modifier only works with vertex-based populations
    VertexBasedCellPopulation<DIM>* p_population =
        dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
    if (p_population == nullptr)
    {
        EXCEPTION("NeighborDependentProliferationModifier requires a VertexBasedCellPopulation.");
    }

    MutableVertexMesh<DIM, DIM>& r_mesh =
        static_cast<MutableVertexMesh<DIM, DIM>&>(p_population->rGetMesh());

    MAKE_PTR(TransitCellProliferativeType, p_transit);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        std::set<unsigned> neighbours = r_mesh.GetNeighbouringElementIndices(elem_index);

        if (neighbours.size() >= mNeighborThreshold)
        {
            cell_iter->SetCellProliferativeType(p_diff);
        }
        else
        {
            cell_iter->SetCellProliferativeType(p_transit);
        }
    }
}

template<unsigned DIM>
void NeighborDependentProliferationModifier<DIM>::SetupSolve(
    AbstractCellPopulation<DIM, DIM>& rCellPopulation,
    std::string outputDirectory)
{
    // Call UpdateAtEndOfTimeStep once at setup so initial state is correct
    UpdateAtEndOfTimeStep(rCellPopulation);
}

template<unsigned DIM>
void NeighborDependentProliferationModifier<DIM>::OutputSimulationModifierParameters(
    out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<NeighborThreshold>" << mNeighborThreshold << "</NeighborThreshold>\n";
    AbstractCellBasedSimulationModifier<DIM, DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class NeighborDependentProliferationModifier<1>;
template class NeighborDependentProliferationModifier<2>;
template class NeighborDependentProliferationModifier<3>;

#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NeighborDependentProliferationModifier)
