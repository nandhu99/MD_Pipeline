import parse_options
import runProcess

all_flags = parse_options.parseConfig2('config.txt')
pdb_file, dssp_file = runProcess.runPreProcess(all_flags)
cg_protein, cg_topol, cg_index, nmap = runProcess.runMaritinize(all_flags, pdb_file, dssp_file)
bool_box, system, system_top, lipids = runProcess.runInsane(all_flags, pdb_file, cg_protein, cg_topol)
system_ndx = runProcess.make_ndx(system, lipids)
topology = runProcess.make_topology(cg_topol, system_top)
em_gro = runProcess.runMinimization(all_flags, system, topology, system_ndx)
equil_gro = runProcess.runEquilibration(all_flags, system, lipids, em_gro, topology, system_ndx)
