import argparse
import math

def distance(p1, p2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))
 
def find_closest_position(positions, target):
    return min(positions, key=lambda pos: distance(pos, target))

def parse_poscar(file_path, verbose=False):
    if verbose:
        print(f"Reading POSCAR file from: {file_path}")
        
    with open(file_path, 'r') as f:
        lines = f.readlines()

    atom_types = lines[5].split()
    atom_counts = list(map(int, lines[6].split()))
    total_atoms = sum(atom_counts)

    if verbose:
        print(f"Parsed atom types: {atom_types}")
        print(f"Parsed atom counts: {atom_counts}")
        print(f"Total number of atoms: {total_atoms}")

    coord_line_index = 7
    while lines[coord_line_index].strip().lower() not in ['direct', 'cartesian']:
        coord_line_index += 1

    coord_type = lines[coord_line_index].strip()

    start_index = coord_line_index + 1
    atomic_positions = [list(map(float, lines[start_index + i].split()[:3])) for i in range(total_atoms)]

    element_positions = {}
    idx = 0
    for atom, count in zip(atom_types, atom_counts):
        element_positions[atom] = atomic_positions[idx:idx + count]
        idx += count

    return {
        "lines": lines,
        "atom_types": atom_types,
        "atom_counts": atom_counts,
        "coord_type": coord_type,
        "element_positions": element_positions
    }

def substitute(parsed_data, substitution, site, target_position=None, verbose=False):
    lines = parsed_data["lines"]
    atom_types = parsed_data["atom_types"].copy()
    atom_counts = parsed_data["atom_counts"].copy()
    coord_type = parsed_data["coord_type"]
    positions = parsed_data["element_positions"]

    if site not in atom_types:
        raise ValueError(f"{site} not found in POSCAR.")

    atom_positions = positions[site]

    if target_position is None:
        target_position = [0.5, 0.5, 0.5]
        if verbose:
            print("No target specified. Using default center: [0.5, 0.5, 0.5]")

    substituted_position = find_closest_position(atom_positions, target_position)
    atom_positions.remove(substituted_position)

    if verbose:
        print(f"Substituting {site} with {substitution}")
        print(f"Selected atom at: {substituted_position} (closest to {target_position})")

    site_index = atom_types.index(site)
    atom_counts[site_index] -= 1

    if atom_counts[site_index] == 0 and verbose:
        print(f"Warning: all {site} atoms have been removed.")

    if substitution in atom_types:
        sub_index = atom_types.index(substitution)
        atom_counts[sub_index] += 1
    elif substitution == 'Va':
        pass
    else:
        atom_types.append(substitution)
        atom_counts.append(1)
        positions[substitution] = []

    if substitution != 'Va':
        positions[substitution].append(substituted_position)
        if verbose:
            print(f"Added {substitution} atom at: {substituted_position}")

    all_positions = []
    for atom in atom_types:
        all_positions.extend(positions.get(atom, []))

    new_lines = lines[:5]
    new_lines.append("  " + "  ".join(atom_types) + "\n")
    new_lines.append("  " + "  ".join(str(n) for n in atom_counts) + "\n")
    new_lines.append(coord_type + "\n")
    for pos in all_positions:
        new_lines.append(f"  {pos[0]:.16f}  {pos[1]:.16f}  {pos[2]:.16f}\n")

    new_lines[0] = f'{substitution}_{site} defect {" ".join(f"{x:.16f}" for x in substituted_position)}\n'
    print(f"({substitution} → {site}, {substituted_position})")

    filename = f"{substitution}_{site}_POSCAR"
    return new_lines, filename

def write_poscar(new_lines, filename, verbose=False):
    with open(filename, 'w') as f:
        f.writelines(new_lines)

    if verbose:
        print(f"Written updated POSCAR to '{filename}'")

def main():
    parser = argparse.ArgumentParser(
        description="Substitute or remove atoms in a VASP POSCAR file."
    )
    parser.add_argument("substitution", type=str, help="Atom to substitute in (or 'Va' for vacancy)")
    parser.add_argument("site", type=str, help="Atom type whose site will be replaced")
    parser.add_argument("-f", "--file", type=str, default="POSCAR", help="Input POSCAR file")
    parser.add_argument("-o", "--output", type=str, default=None, help="Output filename")
    parser.add_argument("--target", nargs=3, type=float, metavar=('X', 'Y', 'Z'),
                        help="Target fractional coordinate to find closest atom to")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()
    target_position = args.target if args.target else None

    parsed_data = parse_poscar(args.file, verbose=args.verbose)
    new_lines, filename = substitute(parsed_data, args.substitution, args.site,
                                     target_position=target_position, verbose=args.verbose)
    output_file = args.output if args.output else filename
    write_poscar(new_lines, output_file, verbose=args.verbose)

if __name__ == "__main__":
    main()
