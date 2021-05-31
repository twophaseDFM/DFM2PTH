// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The spatial parameters for the fracture sub-domain in the exercise
 *        on two-phase flow in fractured porous media.
 */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_SPATIALPARAMS_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_SPATIALPARAMS_HH

#include <dumux/io/grid/griddata.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The spatial params the two-phase facet coupling test
 */
template< class FVGridGeometry, class Scalar >
class FractureSpatialParams
: public FVSpatialParams< FVGridGeometry, Scalar, FractureSpatialParams<FVGridGeometry, Scalar> >
{
    using ThisType = FractureSpatialParams< FVGridGeometry, Scalar >;
    using ParentType = FVSpatialParams< FVGridGeometry, Scalar, ThisType >;

    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

    // we identify those fractures as barriers, that have a domain marker
    // of 2. This is what is set in the grid file (see grids/complex.geo)
    static constexpr int barriersDomainMarker = 2;

public:
    //! export the type used for permeabilities
    using PermeabilityType = Scalar;

    //! the constructor
    FractureSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                          std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                          const std::string& paramGroup)
    : ParentType(fvGridGeometry)
    , gridDataPtr_(gridData)
    , pcKrSwCurve_("Fracture.SpatialParams")
    , barrierPcKrSwCurve_("Fracture.SpatialParams.Barrier")
    {
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
        permeability1_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability1");
        permeability2_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability2");
        permeability3_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability3");
        permeability4_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability4");
        permeability5_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability5");
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    template< class ElementSolution >
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (getElementDomainMarker(element) == 1)
        	return permeability1_;
        else if(getElementDomainMarker(element) == 2)
        	return permeability2_;
        else if(getElementDomainMarker(element) == 3)
        	return permeability3_;
        else if(getElementDomainMarker(element) == 4)
        	return permeability4_;
        else
        	return permeability5_;
    }

    //! Return the porosity
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law for the sub-control volume
     *
     * \param element The current finite element
     * \param scv The sub-control volume
     * \param elemSol The current element solution
     */
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {  
        return makeFluidMatrixInteraction(pcKrSwCurve_); 
    }

    //! Water is the wetting phase
    template< class FluidSystem >
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        // we set water as the wetting phase here
        // which is phase0Idx in the H2oN2 fluid system
        return FluidSystem::phase0Idx;
    }

    //! returns the domain marker for an element
    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

private:
    //! pointer to the grid data (contains domain markers)
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;

    Scalar porosity_;
    Scalar porosityBarrier_;
    PermeabilityType permeability1_,permeability2_,permeability3_,permeability4_,permeability5_ ;
    PermeabilityType permeabilityBarrier_;
    const PcKrSwCurve pcKrSwCurve_;
    const PcKrSwCurve barrierPcKrSwCurve_;
};

} // end namespace Dumux

#endif
