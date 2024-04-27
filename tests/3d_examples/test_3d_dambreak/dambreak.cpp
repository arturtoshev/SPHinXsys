#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// general parameters for geometry
Real dx = 0.025;              // particle spacing
Real BW = dx * 4;             // boundary width
Real DL = 5.366;              // tank length
Real DH = 2.0;                // tank height
Real DW = 2.0;                // tank width
Real LL = 2.0;                // liquid length 0.4  / 2
Real LH = 0.7;                // liquid height 0.3  / 1 
Real LW = 2.0;                // liquid width 0.67  / 0.5
// TODO: DW, LW to 2.0

// for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;

// Shape parameters later overriden by commandline arguments
Real A_top_x = 0.0;
Real A_top_z = 0.0;
Real T_top_x = 0.0;
Real T_top_z = 0.0;
Real S_top_x = 0.0;
Real S_top_z = 0.0;
Real A_front_z = 0.0;
Real A_front_y = 0.0;
Real T_front_z = 0.0;
Real T_front_y = 0.0;
Real S_front_z = 0.0;
Real S_front_y = 0.0;
std::default_random_engine generator;
std::normal_distribution<Real> distribution(0.0, dx / 10.0);

//	define the water block generator
class WaterBlockGenerator : public ParticleGenerator<Base>
{
  public:
    explicit WaterBlockGenerator(SPHBody &sph_body) : ParticleGenerator<Base>(sph_body){};
    virtual void initializeGeometricVariables() override
    {
        // add fluid points
        for (Real x = 0.5 * dx; x <= LL + A_front_y + A_front_z; x += dx)  // length
        {
            for (Real z = 0.5 * dx; z <= LW; z += dx)  // width
            {
                Real water_height = LH + A_top_x * sin(2.0 * M_PI * (T_top_x * x / LL + S_top_x)) 
                                       + A_top_z * sin(2.0 * M_PI * (T_top_z * z / LW + S_top_z));
                for (Real y = 0.5 * dx; y <= LH + A_top_x + A_top_z; y += dx)  // height
                {
                    Real water_length = LL + A_front_z * sin(2.0 * M_PI * (T_front_z * z / LW + S_front_z))
                                           + A_front_y * sin(2.0 * M_PI * (T_front_y * y / LH + S_front_y));
                    if (y <= water_height && x <= water_length) {
                        Real px = x + distribution(generator);
                        Real py = y + distribution(generator);
                        Real pz = z + distribution(generator);
                        initializePositionAndVolumetricMeasure(Vecd(px, py, pz), dx * dx * dx);
                    }
                }
            }
        }
    }
};

//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
        Vecd halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Transform translation_wall(halfsize_inner);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_outer);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_inner);
    }
};

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //  Pass command line arguments to the program
    //  Example run: nohup ./bin/test_3d_dambreak 20.0 0.15 0.15 0.5 1 0 0 ./out_dir_1 42 >> nohup.out 2>&1 &
    // ./bin/test_3d_dambreak 0.05 0.15 0.15 1 0.5 0 0.75 /home/atoshev/code/sphinxsys_vol/out_dir_3 42 0 0.15 1 0.25 0 0.75
    //----------------------------------------------------------------------
    // Duration of the simulation
    Real end_time = atof(av[1]);  // 20.0
    // Shape parameters
    A_top_x = atof(av[2]);  // [0.0, 0.15]  - amplitude of the top wave along x
    A_top_z = atof(av[3]);  // [0.0, 0.15]
    T_top_x = atof(av[4]);  // [0.25, 2]  - period of the top wave along x
    T_top_z = atof(av[5]);  // [0.25, 2]
    S_top_x = atof(av[6]);  // [0.0, 1]  - phase shift of the top wave along x
    S_top_z = atof(av[7]);  // [0.0, 1]
    // Output folder
    std::string output_folder_ = av[8];   // "./output"; - for writing vtk files
    fs::remove_all(output_folder_);
    fs::create_directory(output_folder_);
    // Seed for additive Gaussian noise to initial positions
    unsigned int seed = atoi(av[9]);  // [42[]
    generator.seed(seed);
    // Shape parameters
    A_front_z = atof(av[10]);  // [0.0, 0.15]  - amplitude of the front wave along z
    A_front_y = atof(av[11]);  // [0.0, 0.15]
    T_front_z = atof(av[12]);  // [0.25, 2]  - period of the front wave along z
    T_front_y = atof(av[13]);  // [0.25, 2]
    S_front_z = atof(av[14]);  // [0.0, 1]  - phase shift of the front wave along z
    S_front_y = atof(av[15]);  // [0.0, 1]
    std::cout << std::fixed << std::setprecision(3) << "Shape parameters: "
              << "\n  A_top_x=" << A_top_x << ", A_top_z=" << A_top_z << ", T_top_x=" << T_top_x << ", T_top_z=" << T_top_z << ", S_top_x=" << S_top_x << ", S_top_z=" << S_top_z 
              << "\n  A_front_z=" << A_front_z << ", A_front_y=" << A_front_y << ", T_front_z=" << T_front_z << ", T_front_y=" << T_front_y << ", S_front_z=" << S_front_z << ", S_front_y=" << S_front_y
              << "\n  End time=" << end_time << ", Output folder=" << output_folder_ << ", Seed=" << seed << std::endl;
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DH + BW, DW + BW));
    SPHSystem sph_system(system_domain_bounds, dx);
    int argc = 1; char* argv[] = { av[0] };
    sph_system.handleCommandlineOptions(argc, argv)->setIOEnvironment();
    sph_system.getIOEnvironment().output_folder_ = output_folder_;
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, "WaterBody");
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles(WaterBlockGenerator(water_block));
    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_wall_contact);
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce> constant_gravity(water_block, gravity);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_wall_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_water_block_states(sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;  // sph_system.RestartStep();
    int screen_output_interval = 100;
    int output_interval = 100;  // Real output_interval = end_time / 20.0;
    Real dt = 0.0008; // default acoustic time step sizes
    Real Dt = 0.004;  // get_fluid_advection_time_step_size.exec();
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_water_block_states.writeToFile(0);
    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = " 
            << GlobalStaticVariables::physical_time_ << "	Dt = " << Dt << "	dt = " << dt << "\n";
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < (end_time + 0.00000001))
    {
        update_density_by_summation.exec();

        pressure_relaxation.exec(dt); density_relaxation.exec(dt); GlobalStaticVariables::physical_time_ += dt;
        pressure_relaxation.exec(dt); density_relaxation.exec(dt); GlobalStaticVariables::physical_time_ += dt;
        pressure_relaxation.exec(dt); density_relaxation.exec(dt); GlobalStaticVariables::physical_time_ += dt;
        pressure_relaxation.exec(dt); density_relaxation.exec(dt); GlobalStaticVariables::physical_time_ += dt;
        pressure_relaxation.exec(dt); density_relaxation.exec(dt); GlobalStaticVariables::physical_time_ += dt;

        water_block.updateCellLinkedListWithParticleSort(100);
        water_block_complex.updateConfiguration();
        number_of_iterations = number_of_iterations + 5;


        if (number_of_iterations % screen_output_interval == 0)
        {
            std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = " 
                      << GlobalStaticVariables::physical_time_ << "	Dt = " << Dt << "	dt = " << dt << "\n";
        }

        if (number_of_iterations % output_interval == 0)
        {
            TickCount t2 = TickCount::now();
            write_water_block_states.writeToFile(number_of_iterations);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;  // time for writing files
        }
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval; // total time exluding writing files
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    return 0;
}
