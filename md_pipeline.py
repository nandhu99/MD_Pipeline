import parse_options
import runProcess

all_flags = parse_options.parseConfig2('config.txt')
clean_pdb = runProcess.runPreProcess(all_flags)
Sep_chains, processed_pdb, dssp_file = runProcess.runMemembed(clean_pdb, all_flags)
cg_protein, cg_topol, cg_index, nmap = runProcess.runMaritinize(all_flags, processed_pdb, dssp_file)
multi_prot, total_prot = runProcess.multiplyProtein(all_flags, cg_protein)
bool_box, system, system_top, lipids = runProcess.runInsane(all_flags, clean_pdb, multi_prot)
system_ndx = runProcess.make_ndx(system, lipids)
topology = runProcess.make_topology(clean_pdb, cg_topol, system_top, total_prot)
em_gro, em_tpr = runProcess.runMinimization(all_flags, system, topology, system_ndx)
equil_gro, equil_tpr = runProcess.runEquilibration(all_flags, system, em_gro, topology, system_ndx)
