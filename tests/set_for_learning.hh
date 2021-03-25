//这个文件被用于掌握一个dynamics模块中每个部分的具体含义

/* *************** Class BGKdynamics *********************************************** */

template<typename T, template<typename U> class Descriptor>
//
int BGKdynamics<T,Descriptor>::id =
    meta::registerGeneralDynamics<T,Descriptor,BGKdynamics<T,Descriptor> >("BGK");

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, template<typename U> class Descriptor>
BGKdynamics<T,Descriptor>::BGKdynamics(T omega_ )
    : IsoThermalBulkDynamics<T,Descriptor>(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
BGKdynamics<T,Descriptor>::BGKdynamics(HierarchicUnserializer& unserializer)
    : IsoThermalBulkDynamics<T,Descriptor>(T())
{
    this->unserialize(unserializer);
}

template<typename T, template<typename U> class Descriptor>
BGKdynamics<T,Descriptor>* BGKdynamics<T,Descriptor>::clone() const {
    return new BGKdynamics<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int BGKdynamics<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void BGKdynamics<T,Descriptor>::collide (
        Cell<T,Descriptor>& cell,
        BlockStatistics& statistics )
{
    T rhoBar;//对密度的声明
    Array<T,Descriptor<T>::d> j;//数组容器，维度数量
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);


/**static void get_rhoBar_j(Array<T,Descriptor::q> const& f, T& rhoBar, Array<T,Descriptor::d>& j ) {
    rhoBar = get_rhoBar(f);
    get_j(f, j);
}*/

/**static T get_rhoBar(Array<T,Descriptor::q> const& f) {
    T rhoBar = f[0];
    for (plint iPop=1; iPop < Descriptor::q; ++iPop) {
        rhoBar += f[iPop];
    }
    return rhoBar;
}*/

/**static void get_j(Array<T,Descriptor::q> const& f, Array<T,Descriptor::d>& j ) {
    for (int iD=0; iD < Descriptor::d; ++iD) {
        j[iD] = f[0]*Descriptor::c[0][iD];
    }
    for (plint iPop=1; iPop < Descriptor::q; ++iPop) {
        for (int iD=0; iD < Descriptor::d; ++iD) {
            j[iD] += f[iPop]*Descriptor::c[iPop][iD];
        }
    }
}
static T ‪i-nvRho(T ‪rhoBar) {
         return (T)1 / (‪rhoBar + (T)1);
     }
*/


    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());

/**static T bgk_ma2_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_collision(cell.getRawPopulations(), rhoBar, j, omega);
}*/

/**template<typename T, class Descriptor>
struct dynamicsTemplatesImpl {

static T bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr );
    }
    return jSqr*invRho*invRho;
}*/

    if (cell.takesStatistics()) {
        gatherStatistics(statistics, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
void BGKdynamics<T,Descriptor>::collideExternal (
        Cell<T,Descriptor>& cell, T rhoBar,
        Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat )
{
    T uSqr = dynamicsTemplates<T,Descriptor>::bgk_ma2_collision(cell, rhoBar, j, this->getOmega());
    if (cell.takesStatistics()) {
        gatherStatistics(stat, rhoBar, uSqr);
    }
}

template<typename T, template<typename U> class Descriptor>
T BGKdynamics<T,Descriptor>::computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                                T jSqr, T thetaBar) const
{
    T invRho = Descriptor<T>::invRho(rhoBar);
    return dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

template<typename T, template<typename U> class Descriptor>
void BGKdynamics<T,Descriptor>::decomposeOrder0 (
        Cell<T,Descriptor> const& cell, std::vector<T>& rawData ) const
{
    T rhoBar;
    Array<T,Descriptor<T>::d> j;
    momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    rawData[0] = rhoBar;
    j.to_cArray(&rawData[1]);
    
    Array<T,Descriptor<T>::q> fEq;
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq );

    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        rawData[1+Descriptor<T>::d+iPop] =
            cell[iPop] - fEq[iPop];
    }

    int offset = 1+Descriptor<T>::d+Descriptor<T>::q;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        rawData[offset+iExt] = *cell.getExternal(iExt);
    }
}

template<typename T, template<typename U> class Descriptor>
void BGKdynamics<T,Descriptor>::recomposeOrder0 (
        Cell<T,Descriptor>& cell, std::vector<T> const& rawData ) const
{
    T rhoBar = rawData[0];
    Array<T,Descriptor<T>::d> j;
    j.from_cArray(&rawData[1]);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);

    
    Array<T,Descriptor<T>::q> fEq;
    dynamicsTemplates<T,Descriptor>::bgk_ma2_equilibria( rhoBar, Descriptor<T>::invRho(rhoBar), j, jSqr, fEq );
    
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = fEq[iPop] + rawData[1+Descriptor<T>::d+iPop];
    }

    int offset = 1+Descriptor<T>::d+Descriptor<T>::q;
    for (plint iExt=0; iExt<Descriptor<T>::ExternalField::numScalars; ++iExt) {
        *cell.getExternal(iExt) = rawData[offset+iExt];
    }
}