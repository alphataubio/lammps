
import re

atm_parameters_list = [
    ['r_s', 'valency', 'mass', 'r_vdw', 'epsilon', 'gamma', 'r_pi', 'valency_e'],
    ['alpha', 'gamma_w', 'valency_boc', 'p_ovun5', 'not_used', 'chi', 'eta', 'p_hbond'], 
    ['r_pi_pi', 'p_lp2', 'not_used', 'b_o_131', 'b_o_132', 'b_o_133', 'bcut_acks2', 'not_used'], 
    ['p_ovun2', 'p_val3', 'not_used', 'valency_val', 'p_val5', 'rcore2', 'ecore2', 'acore2']]

bnd_parameters_list = [
    ['De_s','De_p','De_pp','p_be1','p_bo5','v13cor','p_bo6','p_ovun1'], 
    ['p_be2','p_bo3','p_bo4','not_used','p_bo1','p_bo2','ovc','not_used']]

ofd_parameters_list = ['D', 'r_vdW', 'alpha', 'r_s', 'r_p', 'r_pp']

ang_parameters_list = ['theta_00', 'p_val1', 'p_val2', 'p_coa1', 'p_val7', 'p_pen1', 'p_val4']

parameters_list = ['V1', 'V2', 'V3', 'p_tor1', 'p_cot1', 'not_used', 'not_used']

hbd_parameters_list = ['r0_hb', 'p_hb1', 'p_hb2', 'p_hb3']


# reaxff-komissarov2021
 
with open('/Users/mitch/Dropbox/lammps/lammps-alphataubio/potentials/reaxff-zhang2018.ff', 'r') as f:

    # HEADER LINE (AS IS)
    print( f.readline().rstrip('\n') )
 
    # GEN BLOCK
    re_general = r'^(\s*\-?[0-9]+\.[0-9]+)\s*(!?[^!]*$)'
    second_line = f.readline()
    num_general_parameters = int(second_line.split()[0])
    print(f'{num_general_parameters:2}       ! Number of general parameters')

    for i in range(num_general_parameters):
        
        if( not (match := re.search( re_general, f.readline().rstrip('\n') ) )):
            raise Exception(f"Unable to parse GEN block at line {i+3}.")

        #print('{:8.4f} {}'.format(float(match.group(1)), match.group(2)))

    # ATM BLOCK
    re_atm = r'^([A-Za-z]+)(?:\s+\-?[0-9]+\.[0-9]+){32}\s*'
    num_atoms = int(f.readline().split()[0])
    f.readline()
    f.readline()
    f.readline()
    print('{:<3}  ! NUM_ATOMS {}, {}, {}, {}, {}, {}, {}, {}'.format(num_atoms, *atm_parameters_list[0]))
    print('                 {}, {}, {}, {}, {}, {}, {}, {}'.format(*atm_parameters_list[1]))
    print('                 {}, {}, {}, {}, {}, {}, {}, {}'.format(*atm_parameters_list[2]))
    print('                 {}, {}, {}, {}, {}, {}, {}, {}'.format(*atm_parameters_list[3]))

    ff_string = f.read()

    for i in range(num_atoms):
        
        if( not (match := re.match( re_atm, ff_string, flags=re.MULTILINE|re.DOTALL)) ):
            raise Exception("Unable to parse ATM block, missing space between numbers perhaps ?")

        tokens = match.group(0).split()
        tokens_formatted = [f'{match.group(1):2}'] + [' {:8.4f}'.format(float(t)) for t in tokens[1:]]
        tokens_formatted.insert(9, "\n  ")
        tokens_formatted.insert(18, "\n  ")
        tokens_formatted.insert(27, "\n  ")
        #print( ''.join(tokens_formatted))
        ff_string = ff_string[len(match.group(0)):]

    # BND BLOCK
    re_bnd = r'^\s*([0-9]+)\s+([0-9]+)(?:\s+\-?[0-9]+\.[0-9]+){16}\s*'
    num_bonds = int(ff_string[:10].split()[0])
    print('{:<4} ! NUM_BONDS {}, {}, {}, {}, {}, {}, {}, {}'.format(num_bonds, *bnd_parameters_list[0]))
    print('                 {}, {}, {}, {}, {}, {}, {}, {}'.format(*bnd_parameters_list[1]))
    ff_string = ff_string.split('\n', 2)[2]
    
    for i in range(num_bonds):
        
        if( not (match := re.match( re_bnd, ff_string, flags=re.MULTILINE|re.DOTALL)) ):
            raise Exception("Unable to parse BND block, missing space between numbers perhaps ?")

        tokens = match.group(0).split()
        tokens_formatted = [f'{match.group(1):2} {match.group(2):2}'] + [' {:8.4f}'.format(float(t)) for t in tokens[2:]]
        tokens_formatted.insert(9, "\n     ")
        #print( ''.join(tokens_formatted))
        ff_string = ff_string[len(match.group(0)):]

    # OFD BLOCK
    re_ofd = r'^\s*([0-9]+)\s+([0-9]+)(?:\s+\-?[0-9]+\.[0-9]+){6}\s*'
    num_ofd = int(ff_string[:10].split()[0])
    print('{:<4} ! NUM_OFF_DIAGONALS {}, {}, {}, {}, {}, {}'.format(num_ofd, *ofd_parameters_list))
    ff_string = ff_string.split('\n', 1)[1]
    
    for i in range(num_ofd):
        
        if( not (match := re.match( re_ofd, ff_string, flags=re.MULTILINE|re.DOTALL)) ):
            raise Exception("Unable to parse OFD block, missing space between numbers perhaps ?")

        tokens = match.group(0).split()
        tokens_formatted = [f'{match.group(1):2} {match.group(2):2}']
        tokens_formatted += [' {:8.4f}'.format(float(t)) for t in tokens[2:]]
        #print( ''.join(tokens_formatted))
        ff_string = ff_string[len(match.group(0)):]

    # ANG BLOCK
    re_ang = r'^\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)(?:\s+\-?[0-9]+\.[0-9]+){7}\s*'
    num_ang = int(ff_string[:10].split()[0])
    print('{:<4} ! NUM_ANGLES {}, {}, {}, {}, {}, {}, {}'.format(num_ang, *ang_parameters_list))
    ff_string = ff_string.split('\n', 1)[1]
    
    for i in range(num_ang):
        
        if( not (match := re.match( re_ang, ff_string, flags=re.MULTILINE|re.DOTALL)) ):
            raise Exception("Unable to parse ANG block, missing space between numbers perhaps ?")

        #print( match.group(0) )
        tokens = match.group(0).split()
        tokens_formatted = [f'{match.group(1):2} {match.group(2):2} {match.group(3):2}'] 
        tokens_formatted += [' {:8.4f}'.format(float(t)) for t in tokens[3:]]
        print( ''.join(tokens_formatted))
        ff_string = ff_string[len(match.group(0)):]


