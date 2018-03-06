# Modified based on code from http://www.pymolwiki.org/index.php/FindSurfaceResidues
# Author: D. Allan Drummond and Nicholas W. VanKuren
import pymol, sys, random
from pymol import cmd

# From http://www.pymolwiki.org/index.php/Aa_codes
# one_letter["SER"] will now return "S"
one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', \
'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',    \
'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',    \
'GLY':'G', 'PRO':'P', 'CYS':'C'}

# From Tien et al. 2012 (unpublished), http://arxiv.org/pdf/1211.4251v2.pdf
# Maximum theoretical surface area of residue X, in Gly-X-Gly tripeptide, in square angstroms
maximum_surface_area = {
'ALA':129,'CYS':158,'ASP':193,'GLU':223,
'PHE':239,'GLY':104,'HIS':209,'ILE':197,
'LYS':237,'LEU':201,'MET':218,'ASN':195,
'PRO':159,'GLN':224,'ARG':274,'SER':151,
'THR':172,'VAL':174,'TRP':282,'TYR':263}

def getRelativeSolventAccessibility(objSel="(all)", output_file="/home/justin/Desktop/graduateschoolclasses/Winter_2018/Evolution_of_Biological_Molecules/HW5/solvent.txt", dotdensity=2, solventradius=1.4):
    """
    getSolventAccessibleSurfaceArea
        For each residue in selection, calculate solvent-accessible surface area (SASA) 
        and the relative solvent accessibility (RSA), which is the SASA divided by the 
        maximum achievable SASA
 
    PARAMS
        objSel (string)
            the object or selection in which to calculate SASA
            DEFAULT: (all)
            
    RETURNS
        A list of (resi, chain, resn, SASA, maximum SASA, and RSA
     """
     # Set up output.
    if output_file is None:
         output = sys.stdout
    else:
        output = open(os.path.expanduser(output_file),'w')
     
    # Make a temporary object name
    tmpObjName = "tmp" + str(random.randint(0,1000000))
    # Create a new selection that we can treat as we please; include all atoms by specifying "and polymer"
    cmd.create(tmpObjName, objSel) # + " and polymer")
    # Add hydrogens
    cmd.h_add(tmpObjName)
    
    # Save current values of environmental variables, so that we can set them during this routine.
    current_dot_solvent = cmd.get("dot_solvent")
    current_dot_density = cmd.get("dot_density")
    current_solvent_radius = cmd.get("solvent_radius")
    
    # Set values for assessing surface area
    cmd.set("dot_solvent", 1) # solvent-exposed area
    cmd.set("dot_density", dotdensity)
    cmd.set("solvent_radius", solventradius) # radius, in angstroms, of solvent molecule
    
    # Compute the solvent-exposed area for each residue and load the result into the residue's b-factor
    cmd.get_area(tmpObjName, load_b=True)
    # Compile the b-factors into residue-specific values
    resinfo = {}
    resinfo["sasa_list"]=[]
    cmd.iterate(tmpObjName, 'sasa_list.append((chain,resi,resn,b))', space=resinfo)
    
    solvent_asa_dict = {}
    # Read out the surface areas, per atom, and sum them into surface areas per residue
    for (chain,resi,resn,sasa) in resinfo["sasa_list"]:
        if maximum_surface_area.has_key(resn): # If this is a recognized residue
            cur_sasa = solvent_asa_dict.get((chain,resi,resn), 0.0) # get existing surface area for this residue, or 0.0
            solvent_asa_dict[(chain,resi,resn)] = cur_sasa + sasa # add to residue's surface area.
    
    # Write out the results, sorted by residue number
    results = []
    output.write("site\tchain\taa\tsasa\tmax.sasa\trsa\n")
    for (chain,resi,resn) in sorted(solvent_asa_dict.keys(), key=lambda x: int(x[1])):
        sasa = solvent_asa_dict[(chain,resi,resn)]
        prop_solv_exposed = sasa/maximum_surface_area[resn]
        output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (resi, chain, one_letter[resn], sasa, maximum_surface_area[resn], prop_solv_exposed))
        results.append((resi, chain, resn, sasa, maximum_surface_area[resn], prop_solv_exposed))
    
    # Restore values we changed
    cmd.set("dot_solvent", current_dot_solvent)
    cmd.set("dot_density", current_dot_density)
    cmd.set("solvent_radius", current_solvent_radius)
    # Remove the temporary object
    cmd.delete(tmpObjName)

    # Return results
    return results

cmd.extend("getRelativeSolventAccessibility", getRelativeSolventAccessibility)

cmd.do('''getRelativeSolventAccessibility ''')
