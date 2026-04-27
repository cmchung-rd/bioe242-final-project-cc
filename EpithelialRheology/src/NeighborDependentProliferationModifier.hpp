#ifndef NEIGHBORDEPENDENTPROLIFERATIONMODIFIER_HPP_
#define NEIGHBORDEPENDENTPROLIFERATIONMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A simulation modifier that switches cell proliferative types based on
 * the number of neighbouring elements in a vertex mesh.
 *
 * Cells with >= mNeighborThreshold neighbours are set to
 * DifferentiatedCellProliferativeType (arrest), while cells below the
 * threshold are set to TransitCellProliferativeType (proliferating).
 *
 * This confines division to frontier cells adjacent to free space,
 * preventing interior overcrowding.
 */
template<unsigned DIM>
class NeighborDependentProliferationModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:

    /** Cells with >= this many neighbours are arrested. */
    unsigned mNeighborThreshold;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM> >(*this);
        archive & mNeighborThreshold;
    }

public:

    /**
     * Constructor.
     */
    NeighborDependentProliferationModifier();

    /**
     * Destructor.
     */
    virtual ~NeighborDependentProliferationModifier() = default;

    /**
     * Set the neighbour threshold. Cells with >= this many neighbours are arrested.
     * @param threshold the neighbour count threshold
     */
    void SetNeighborThreshold(unsigned threshold);

    /**
     * @return the current neighbour threshold
     */
    unsigned GetNeighborThreshold() const;

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     * Iterates over all cells and sets their proliferative type based on
     * how many neighbouring elements they have.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * @param rParamsFile the file stream to write parameters to
     */
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NeighborDependentProliferationModifier)

#endif // NEIGHBORDEPENDENTPROLIFERATIONMODIFIER_HPP_
