'''
Created on July 9, 2016
@author: Andrew Abi-Mansour
'''

#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# -------------------------------------------------------------------------

import os, sys, numpy, site
import subprocess
from setuptools import setup, find_packages
from Cython.Build import cythonize
from distutils.command.install import install

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def link(dir):
	txt="""[Desktop Entry]
Version=1.0
Type=Application
Name=PyGran
Comment=graphical User Interface for PyGran
Exec=python -m PyGran
Icon={}/Gran.png
Terminal=false
MimeType=image/x-foo
NotShowIn=KDE""".format(dir)

	with open('PyGran.desktop', 'w') as fp:
		fp.write(txt)

	os.system('chmod +x PyGran.desktop')
	os.system('mv PyGran.desktop ~/.local/share/applications')

class LIGGGHTS(install):


    def execute(self, cmd, cwd):
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, cwd=cwd, shell=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line 
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    def print_progress(self, iteration, prefix='', suffix='', decimals=1, total = 100):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            bar_length  - Optional  : character length of bar (Int)
        """
        str_format = "{0:." + str(decimals) + "f}"
        percents = str_format.format(100 * (iteration / float(total)))
        sys.stdout.write('\r %s%s %s' % (percents, '%', suffix))
        sys.stdout.flush()

    def do_pre_install_stuff(self):
        os.chdir('../liggghts-dem/src')
        os.system('make clean-all')
        os.system('make no-all')
        os.system('make yes-rigid')
        os.system('make makeshlib')
        
        files = 'angle.cpp angle_hybrid.cpp atom.cpp atom_map.cpp atom_vec_atomic.cpp atom_vec_charge.cpp atom_vec.cpp atom_vec_ellipsoid.cpp atom_vec_hybrid.cpp atom_vec_line.cpp atom_vec_sph.cpp atom_vec_sphere.cpp atom_vec_sphere_w.cpp atom_vec_sph_var.cpp bond.cpp bond_hybrid.cpp bounding_box.cpp cfd_datacoupling.cpp cfd_datacoupling_file.cpp cfd_datacoupling_mpi.cpp cfd_regionmodel_none.cpp change_box.cpp citeme.cpp comm.cpp compute_atom_molecule.cpp compute_bond_local.cpp compute_centro_atom.cpp compute_cluster_atom.cpp compute_cna_atom.cpp compute_com.cpp compute_com_molecule.cpp compute_contact_atom.cpp compute_coord_atom.cpp compute.cpp compute_displace_atom.cpp compute_erotate_multisphere.cpp compute_erotate_sphere_atom.cpp compute_erotate_sphere.cpp compute_group_group.cpp compute_gyration.cpp compute_gyration_molecule.cpp compute_inertia_molecule.cpp compute_ke_atom.cpp compute_ke.cpp compute_ke_multisphere.cpp compute_msd.cpp compute_msd_molecule.cpp compute_nparticles_tracer_region.cpp compute_pair_gran_local.cpp compute_pe_atom.cpp compute_pe.cpp compute_pressure.cpp compute_property_atom.cpp compute_property_local.cpp compute_property_molecule.cpp compute_rdf.cpp compute_reduce.cpp compute_reduce_region.cpp compute_rigid.cpp compute_slice.cpp compute_stress_atom.cpp compute_surface.cpp compute_temp.cpp contact_models.cpp container_base.cpp create_atoms.cpp create_box.cpp custom_value_tracker.cpp delete_atoms.cpp delete_bonds.cpp dihedral.cpp dihedral_hybrid.cpp displace_atoms.cpp domain.cpp dump_atom_vtk.cpp dump.cpp dump_custom.cpp dump_custom_vtk.cpp dump_decomposition_vtk.cpp dump_euler_vtk.cpp dump_image.cpp dump_local.cpp dump_local_gran_vtk.cpp dump_mesh_stl.cpp dump_mesh_vtk.cpp dump_movie.cpp dump_xyz.cpp error.cpp finish.cpp fix_adapt.cpp fix_addforce.cpp fix_ave_atom.cpp fix_ave_correlate.cpp fix_ave_euler.cpp fix_aveforce.cpp fix_ave_histo.cpp fix_ave_spatial.cpp fix_ave_time.cpp fix_base_liggghts.cpp fix_box_relax.cpp fix_buoyancy.cpp fix_cfd_coupling_convection.cpp fix_cfd_coupling_convection_impl.cpp fix_cfd_coupling_convection_species.cpp fix_cfd_coupling.cpp fix_cfd_coupling_force.cpp fix_cfd_coupling_force_implicit.cpp fix_check_timestep_gran.cpp fix_check_timestep_sph.cpp fix_contact_history.cpp fix_contact_history_mesh.cpp fix.cpp fix_deform_check.cpp fix_deform.cpp fix_diam_max.cpp fix_drag.cpp fix_dt_reset.cpp fix_efield.cpp fix_enforce2d.cpp fix_external.cpp fix_fiber_spring_simple.cpp fix_freeze.cpp fix_gravity.cpp fix_heat_gran_conduction.cpp fix_heat_gran.cpp fix_insert.cpp fix_insert_pack.cpp fix_insert_rate_region.cpp fix_insert_stream.cpp fix_lb_coupling_onetoone.cpp fix_lineforce.cpp fix_massflow_mesh.cpp fix_mesh.cpp fix_mesh_surface.cpp fix_mesh_surface_stress.cpp fix_mesh_surface_stress_servo.cpp fix_minimize.cpp fix_momentum.cpp fix_move.cpp fix_move_mesh.cpp fix_multisphere.cpp fix_neighlist_mesh.cpp fix_nve.cpp fix_nve_limit.cpp fix_nve_noforce.cpp fix_nve_sph.cpp fix_nve_sphere.cpp fix_nve_sph_stationary.cpp fix_particledistribution_discrete.cpp fix_planeforce.cpp fix_print.cpp fix_property_atom.cpp fix_property_atom_timetracer.cpp fix_property_atom_tracer.cpp fix_property_atom_tracer_stream.cpp fix_property_global.cpp fix_read_restart.cpp fix_region_variable.cpp fix_respa.cpp fix_scalar_transport_equation.cpp fix_setforce.cpp fix_sph.cpp fix_sph_density_continuity.cpp fix_sph_density_corr.cpp fix_sph_density_summation.cpp fix_sph_pressure.cpp fix_spring.cpp fix_spring_rg.cpp fix_spring_self.cpp fix_store.cpp fix_store_force.cpp fix_store_state.cpp fix_template_multiplespheres.cpp fix_template_multisphere.cpp fix_template_sphere.cpp fix_viscous.cpp fix_wall.cpp fix_wall_gran.cpp fix_wall_region.cpp fix_wall_region_sph.cpp fix_wall_sph.cpp force.cpp global_properties.cpp granular_pair_style.cpp granular_styles.cpp granular_wall.cpp group.cpp image.cpp improper.cpp improper_hybrid.cpp input.cpp input_mesh_tet.cpp input_mesh_tri.cpp input_multisphere.cpp integrate.cpp irregular.cpp kspace.cpp lammps.cpp lattice.cpp library_cfd_coupling.cpp library.cpp  math_extra.cpp memory.cpp mesh_mover.cpp min_cg.cpp min.cpp minimize.cpp min_linesearch.cpp modified_andrew.cpp modify.cpp modify_liggghts.cpp multisphere.cpp neigh_bond.cpp neighbor.cpp neigh_derive.cpp neigh_dummy.cpp neigh_full.cpp neigh_gran.cpp neigh_half_bin.cpp neigh_half_multi.cpp neigh_half_nsq.cpp neigh_list.cpp neigh_request.cpp neigh_respa.cpp neigh_stencil.cpp output.cpp pair.cpp pair_gran.cpp pair_gran_proxy.cpp pair_hybrid.cpp pair_hybrid_overlay.cpp pair_soft.cpp pair_sph_artvisc_tenscorr.cpp pair_sph.cpp particleToInsert.cpp particleToInsert_multisphere.cpp procmap.cpp properties.cpp property_registry.cpp random_mars.cpp random_park.cpp read_data.cpp read_dump.cpp reader.cpp reader_native.cpp reader_xyz.cpp read_restart.cpp region_block.cpp region_cone.cpp region.cpp region_cylinder.cpp region_intersect.cpp region_mesh_tet.cpp region_neighbor_list.cpp region_plane.cpp region_prism.cpp region_sphere.cpp region_union.cpp region_wedge.cpp replicate.cpp respa.cpp rotate.cpp run.cpp set.cpp special.cpp tet_mesh.cpp thermo.cpp timer.cpp tri_mesh.cpp tri_mesh_planar.cpp universe.cpp update.cpp variable.cpp velocity.cpp verlet.cpp write_data.cpp write_dump.cpp write_restart.cpp'
        files = files.split()

        count = 0
        print('Compiling LIGGGHTS as a shared library')
        for path in self.execute(cmd='make -f Makefile.shlib openmpi', cwd='/home/levnon/Desktop/liggghts-dem/src'):
            count +=1
            self.print_progress(count, prefix = 'Progress:', suffix = 'Complete', total = len(files) * 2.025)

        sys.stdout.write('\nInstallation complete\n')

    def run(self):
        self.do_pre_install_stuff()
        
setup(
    name = "PyGran",
    version = "0.0.1",
    author = "Andrew Abi-Mansour",
    author_email = "andrew.gaam@gmail.com",
    description = ("A toolbox for rapid quantitative analysis of granular systems"),
    license = "GNU v3",
    keywords = "Discrete Element Method, Granular Materials",
    url = "https://github.com/Andrew-AbiMansour/PyGran",
    packages=find_packages(),
    package_dir={'PyGran': 'PyGran'},
    package_data={'': ['.config', 'GUI/Icons/*.png']},
    install_requires=['numpy', 'pyvtk', 'pytools>=2011.2', 'pygmsh'],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 1 - Alpha",
        "Topic :: Utilities",
        "License :: GNU License",
    ],
    cmdclass={'build_liggghts': LIGGGHTS},
    zip_safe=False,
    ext_modules=cythonize("PyGran/Analyzer/*.pyx"),
	include_dirs=[numpy.get_include()]
)

if sys.argv[1] == 'install':

    sys.path.remove(os.getcwd()) # remove current path to make sure PyGran is imported from elsewhere
    print('Verifying installation ...')
    print('...........................................................................')
    print('...........................................................................')
    print('...........................................................................')
    
    try:
        import PyGran
        link(PyGran.__path__[0] + '/GUI/Icons')
        print('PyGran successfully installed')
    except:
        print('PyGran installation failed ...')
        raise
