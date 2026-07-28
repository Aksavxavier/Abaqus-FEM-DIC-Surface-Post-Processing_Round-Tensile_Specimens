# -*- coding: utf-8 -*-
"""
Export OUTER_CYL_SURF from RB_L0_25_Y_MBW.odb to CSV, VTU and PVD.

Run:
    abaqus python export_RB_L0_25_Y_OUTER_CYL_SURF.py

Model:
    specimen axis = global Y
    axial displacement = U2
    axial logarithmic strain = LE22
    surface = OUTER_CYL_SURF
    element type = C3D8R

Method:
    LE is read at integration points. Each surface-element value is mapped only
    to the nodes of its selected outer face. Contributions from neighboring
    surface elements are averaged at shared surface nodes.
"""

from __future__ import print_function

from odbAccess import openOdb
from abaqusConstants import INTEGRATION_POINT
from collections import OrderedDict
import math
import os


# ============================================================
# USER SETTINGS
# ============================================================

SCRIPT_FOLDER = os.path.dirname(os.path.abspath(__file__))

ODB_PATH = os.path.join(
    SCRIPT_FOLDER,
    "RB_L0_25_Y_MBW.odb"
)

STEP_NAME = "Tensile_Loading"

SURFACE_NAME = "OUTER_CYL_SURF"


# True:
#     Export a sequence of frames.
#
# False:
#     Export only the last frame or REQUESTED_FRAME_ID.
EXPORT_ALL_FRAMES = True


# First frame to export.
FRAME_START = 0


# None:
#     Export until the final available frame.
#
# For a short test:
#     FRAME_END = 2
#
# For all frames:
#     FRAME_END = None
FRAME_END = None


# These settings are used only when:
#
#     EXPORT_ALL_FRAMES = False
#
USE_LAST_FRAME = True

REQUESTED_FRAME_ID = 100


# Export one CSV file for each frame.
EXPORT_CSV = True


OUTPUT_FOLDER = os.path.join(
    SCRIPT_FOLDER,
    "RB_L0_25_Y_OUTER_CYL_SURF_export"
)


PVD_PATH = os.path.join(
    OUTPUT_FOLDER,
    "RB_L0_25_Y_OUTER_CYL_SURF_all_frames.pvd"
)


# ============================================================
# FRAME SELECTION
# ============================================================

def get_frames_to_export(odb):
    """
    Return a list of:

        (frame_id, frame)

    according to the user settings.
    """

    if STEP_NAME not in odb.steps:

        raise RuntimeError(
            "Step not found: %s. Available steps: %s"
            % (
                STEP_NAME,
                list(odb.steps.keys())
            )
        )

    step = odb.steps[STEP_NAME]

    if len(step.frames) == 0:

        raise RuntimeError(
            "Step '%s' contains no frames."
            % STEP_NAME
        )

    if EXPORT_ALL_FRAMES:

        start_id = FRAME_START

        if FRAME_END is None:

            end_id = len(step.frames) - 1

        else:

            end_id = min(
                FRAME_END,
                len(step.frames) - 1
            )

        if start_id < 0 or start_id > end_id:

            raise RuntimeError(
                "Invalid frame range %s to %s. "
                "Available range is 0 to %d."
                % (
                    FRAME_START,
                    FRAME_END,
                    len(step.frames) - 1
                )
            )

        frames = []

        for frame_id in range(
                start_id,
                end_id + 1):

            frames.append(
                (
                    frame_id,
                    step.frames[frame_id]
                )
            )

        return frames

    if USE_LAST_FRAME:

        frame_id = len(step.frames) - 1

    else:

        frame_id = min(
            REQUESTED_FRAME_ID,
            len(step.frames) - 1
        )

    return [
        (
            frame_id,
            step.frames[frame_id]
        )
    ]


# ============================================================
# SURFACE SEARCH
# ============================================================

def surface_name_matches(
        existing_name,
        requested_name):
    """
    Check whether an ODB surface name matches the requested name.
    """

    existing_upper = existing_name.upper()

    requested_upper = requested_name.upper()

    return (
        existing_upper == requested_upper
        or
        existing_upper.endswith(
            "_" + requested_upper
        )
    )


def find_surface_and_instance(assembly):
    """
    Search for OUTER_CYL_SURF first at assembly level and then
    inside every instance.

    Returns:

        surface
        instance_name
    """

    # --------------------------------------------------------
    # Search assembly-level surfaces
    # --------------------------------------------------------

    for key in assembly.surfaces.keys():

        if surface_name_matches(
                key,
                SURFACE_NAME):

            surface = assembly.surfaces[key]

            if len(surface.instanceNames) != 1:

                raise RuntimeError(
                    "Assembly surface must belong "
                    "to exactly one instance."
                )

            instance_name = \
                surface.instanceNames[0]

            return surface, instance_name

    # --------------------------------------------------------
    # Search instance-level surfaces
    # --------------------------------------------------------

    for instance_name, instance in \
            assembly.instances.items():

        for key in instance.surfaces.keys():

            if surface_name_matches(
                    key,
                    SURFACE_NAME):

                surface = \
                    instance.surfaces[key]

                return surface, instance_name

    # --------------------------------------------------------
    # Surface was not found
    # --------------------------------------------------------

    available = []

    for key in assembly.surfaces.keys():

        available.append(
            "ASSEMBLY:" + key
        )

    for instance_name, instance in \
            assembly.instances.items():

        for key in instance.surfaces.keys():

            available.append(
                instance_name +
                ":" +
                key
            )

    raise RuntimeError(
        "Surface '%s' not found. "
        "Available surfaces: %s"
        % (
            SURFACE_NAME,
            available
        )
    )


# ============================================================
# SURFACE ELEMENT-FACE PAIRS
# ============================================================

def iter_surface_element_face_pairs(surface):
    """
    Yield:

        element
        face_name

    for every selected face belonging to the ODB surface.
    """

    elements = surface.elements

    faces = surface.faces

    if len(elements) == 0:

        return

    first = elements[0]

    # --------------------------------------------------------
    # Simple element array
    # --------------------------------------------------------

    if hasattr(
            first,
            "connectivity"):

        for element, face in zip(
                elements,
                faces):

            yield element, str(face)

    # --------------------------------------------------------
    # Grouped element arrays
    # --------------------------------------------------------

    else:

        for element_array, face_array in zip(
                elements,
                faces):

            for element, face in zip(
                    element_array,
                    face_array):

                yield element, str(face)


# ============================================================
# C3D8R FACE CONNECTIVITY
# ============================================================

def get_face_node_labels(
        element,
        face_name):
    """
    Return the node labels belonging to one C3D8/C3D8R face.

    Abaqus element connectivity uses eight nodes.

    Python indices are zero-based.
    """

    connectivity = list(
        element.connectivity
    )

    face_map = {

        "FACE1": [
            0,
            1,
            2,
            3
        ],

        "FACE2": [
            4,
            5,
            6,
            7
        ],

        "FACE3": [
            0,
            4,
            5,
            1
        ],

        "FACE4": [
            1,
            5,
            6,
            2
        ],

        "FACE5": [
            2,
            6,
            7,
            3
        ],

        "FACE6": [
            3,
            7,
            4,
            0
        ]
    }

    if len(connectivity) != 8:

        return []

    if face_name not in face_map:

        return []

    labels = []

    for index in face_map[face_name]:

        labels.append(
            connectivity[index]
        )

    return labels


# ============================================================
# BUILD OUTER SURFACE MESH
# ============================================================

def build_surface_mesh(
        surface,
        instance):
    """
    Build:

        surface node labels
        surface reference coordinates
        selected face connectivity
        surface element-label set
    """

    node_dictionary = {}

    for node in instance.nodes:

        node_dictionary[node.label] = node

    surface_node_labels = set()

    surface_element_labels = set()

    raw_cells = []

    # --------------------------------------------------------
    # Read all selected element faces
    # --------------------------------------------------------

    for element, face_name in \
            iter_surface_element_face_pairs(
                surface):

        face_nodes = get_face_node_labels(
            element,
            face_name
        )

        if len(face_nodes) < 3:

            continue

        surface_element_labels.add(
            element.label
        )

        raw_cells.append(
            face_nodes
        )

        for node_label in face_nodes:

            surface_node_labels.add(
                node_label
            )

    node_labels = sorted(
        surface_node_labels
    )

    if len(node_labels) == 0:

        raise RuntimeError(
            "No nodes found on surface %s."
            % SURFACE_NAME
        )

    # --------------------------------------------------------
    # Build point list
    # --------------------------------------------------------

    node_to_index = {}

    points = []

    for index, node_label in enumerate(
            node_labels):

        if node_label not in node_dictionary:

            raise RuntimeError(
                "Node %d not found in instance %s."
                % (
                    node_label,
                    instance.name
                )
            )

        node_to_index[node_label] = index

        coordinates = \
            node_dictionary[
                node_label
            ].coordinates

        points.append(
            (
                float(coordinates[0]),
                float(coordinates[1]),
                float(coordinates[2])
            )
        )

    # --------------------------------------------------------
    # Convert Abaqus node labels to VTU point indices
    # --------------------------------------------------------

    cells = []

    for face_nodes in raw_cells:

        cell = []

        for node_label in face_nodes:

            cell.append(
                node_to_index[
                    node_label
                ]
            )

        cells.append(cell)

    if len(cells) == 0:

        raise RuntimeError(
            "No cells built for surface %s."
            % SURFACE_NAME
        )

    return (
        node_labels,
        points,
        cells,
        surface_element_labels
    )


# ============================================================
# DISPLACEMENT EXTRACTION
# ============================================================

def make_three_component_vector(data):
    """
    Convert Abaqus vector data to a three-component tuple.
    """

    try:

        number_of_components = len(data)

    except TypeError:

        return (
            0.0,
            0.0,
            0.0
        )

    if number_of_components >= 3:

        return (
            float(data[0]),
            float(data[1]),
            float(data[2])
        )

    if number_of_components == 2:

        return (
            float(data[0]),
            float(data[1]),
            0.0
        )

    if number_of_components == 1:

        return (
            float(data[0]),
            0.0,
            0.0
        )

    return (
        0.0,
        0.0,
        0.0
    )


def get_displacement_dictionary(
        frame,
        instance_name):
    """
    Return:

        node_label -> (U1, U2, U3)
    """

    if "U" not in frame.fieldOutputs:

        raise RuntimeError(
            "U field output is missing."
        )

    displacement_dictionary = {}

    for value in \
            frame.fieldOutputs["U"].values:

        if value.instance is None:

            continue

        if value.instance.name != \
                instance_name:

            continue

        displacement_dictionary[
            value.nodeLabel
        ] = make_three_component_vector(
            value.data
        )

    return displacement_dictionary


# ============================================================
# STRAIN EXTRACTION AND SURFACE-NODE AVERAGING
# ============================================================

def average_le_at_surface_nodes(
        frame,
        surface,
        instance_name):
    """
    Read LE at integration points and average it to the nodes
    of the selected OUTER_CYL_SURF faces.

    Abaqus 3D symmetric tensor ordering:

        index 0 = LE11
        index 1 = LE22
        index 2 = LE33
        index 3 = LE12
        index 4 = LE13
        index 5 = LE23
    """

    if "LE" not in frame.fieldOutputs:

        raise RuntimeError(
            "LE field output is missing."
        )

    le_ip = \
        frame.fieldOutputs["LE"].getSubset(
            position=INTEGRATION_POINT
        )

    # ========================================================
    # First average integration-point values by element
    # ========================================================

    element_accumulator = {}

    for value in le_ip.values:

        if value.instance is None:

            continue

        if value.instance.name != \
                instance_name:

            continue

        if len(value.data) < 6:

            continue

        element_label = value.elementLabel

        if element_label not in \
                element_accumulator:

            element_accumulator[
                element_label
            ] = [
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0
            ]

        for component_index in range(6):

            element_accumulator[
                element_label
            ][component_index] += \
                float(
                    value.data[
                        component_index
                    ]
                )

        element_accumulator[
            element_label
        ][6] += 1

    element_average = {}

    for element_label, values in \
            element_accumulator.items():

        count = values[6]

        if count <= 0:

            continue

        element_average[
            element_label
        ] = (
            values[0] / float(count),
            values[1] / float(count),
            values[2] / float(count),
            values[3] / float(count),
            values[4] / float(count),
            values[5] / float(count)
        )

    # ========================================================
    # Transfer each surface-element value only to the nodes
    # of its selected OUTER_CYL_SURF face
    # ========================================================

    node_accumulator = {}

    for element, face_name in \
            iter_surface_element_face_pairs(
                surface):

        if element.label not in \
                element_average:

            continue

        element_le = \
            element_average[
                element.label
            ]

        face_nodes = get_face_node_labels(
            element,
            face_name
        )

        for node_label in face_nodes:

            if node_label not in \
                    node_accumulator:

                node_accumulator[
                    node_label
                ] = [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0
                ]

            for component_index in range(6):

                node_accumulator[
                    node_label
                ][component_index] += \
                    element_le[
                        component_index
                    ]

            node_accumulator[
                node_label
            ][6] += 1

    # ========================================================
    # Average neighboring surface-element values at each node
    # ========================================================

    nodal_average = {}

    for node_label, values in \
            node_accumulator.items():

        count = values[6]

        if count <= 0:

            continue

        nodal_average[
            node_label
        ] = (
            values[0] / float(count),
            values[1] / float(count),
            values[2] / float(count),
            values[3] / float(count),
            values[4] / float(count),
            values[5] / float(count)
        )

    if len(nodal_average) == 0:

        raise RuntimeError(
            "No LE values could be mapped "
            "to surface %s."
            % SURFACE_NAME
        )

    return nodal_average


# ============================================================
# CYLINDRICAL TRANSFORMATION
# CYLINDER AXIS = GLOBAL Y
# ============================================================

def cylindrical_values(
        reference_point,
        displacement,
        global_le):
    """
    Transform global displacement and strain to cylindrical components.

    Cylinder axis:

        global Y

    Local unit vectors:

        axial:
            e_y = (0, 1, 0)

        radial:
            e_r = (x/r, 0, z/r)

        circumferential:
            e_theta = (-z/r, 0, x/r)

    The reference surface coordinates are used to define the
    radial and circumferential directions.
    """

    x = float(
        reference_point[0]
    )

    z = float(
        reference_point[2]
    )

    radius = math.sqrt(
        x * x +
        z * z
    )

    if radius <= 1.0e-12:

        raise RuntimeError(
            "Cylindrical basis is undefined "
            "at radius zero."
        )

    # --------------------------------------------------------
    # Radial unit vector
    # --------------------------------------------------------

    er_x = x / radius

    er_z = z / radius

    # --------------------------------------------------------
    # Circumferential unit vector
    # --------------------------------------------------------

    et_x = -z / radius

    et_z = x / radius

    # --------------------------------------------------------
    # Global displacement components
    # --------------------------------------------------------

    u1 = displacement[0]

    u2 = displacement[1]

    u3 = displacement[2]

    # --------------------------------------------------------
    # Global logarithmic-strain components
    # --------------------------------------------------------

    le11 = global_le[0]

    le22 = global_le[1]

    le33 = global_le[2]

    le12 = global_le[3]

    le13 = global_le[4]

    le23 = global_le[5]

    # --------------------------------------------------------
    # Local displacement components
    # --------------------------------------------------------

    u_axial = u2

    u_radial = (
        er_x * u1 +
        er_z * u3
    )

    u_circumferential = (
        et_x * u1 +
        et_z * u3
    )

    # --------------------------------------------------------
    # Axial strain
    # --------------------------------------------------------

    le_axial = le22

    # --------------------------------------------------------
    # Radial strain
    # --------------------------------------------------------

    le_radial = (
        er_x * er_x * le11
        +
        er_z * er_z * le33
        +
        2.0 * er_x * er_z * le13
    )

    # --------------------------------------------------------
    # Circumferential strain
    # --------------------------------------------------------

    le_circumferential = (
        et_x * et_x * le11
        +
        et_z * et_z * le33
        +
        2.0 * et_x * et_z * le13
    )

    # --------------------------------------------------------
    # Axial-circumferential tensor shear strain
    # --------------------------------------------------------

    le_ytheta_tensor = (
        et_x * le12
        +
        et_z * le23
    )

    # --------------------------------------------------------
    # Engineering shear strain
    # --------------------------------------------------------

    gamma_ytheta_engineering = (
        2.0 *
        le_ytheta_tensor
    )

    # --------------------------------------------------------
    # Angular position around global Y-axis
    # --------------------------------------------------------

    theta_radians = math.atan2(
        z,
        x
    )

    theta_degrees = (
        theta_radians *
        180.0 /
        math.pi
    )

    return {

        "RADIUS":
            radius,

        "THETA_RAD":
            theta_radians,

        "THETA_DEG":
            theta_degrees,

        "U_AXIAL_Y":
            u_axial,

        "U_RADIAL":
            u_radial,

        "U_CIRCUMFERENTIAL":
            u_circumferential,

        "LE_AXIAL_Y":
            le_axial,

        "LE_RADIAL":
            le_radial,

        "LE_CIRCUMFERENTIAL":
            le_circumferential,

        "LE_YTHETA_TENSOR":
            le_ytheta_tensor,

        "GAMMA_YTHETA_ENGINEERING":
            gamma_ytheta_engineering
    }


# ============================================================
# CSV WRITER
# ============================================================

def write_csv(
        path,
        frame_id,
        frame_value,
        instance_name,
        node_labels,
        points,
        displacements,
        global_le_values,
        local_values):
    """
    Write one row for every outer-surface node.
    """

    header = (
        "Step,"
        "FrameID,"
        "FrameValue,"
        "Instance,"
        "NodeLabel,"
        "X,"
        "Y,"
        "Z,"
        "Radius,"
        "Theta_rad,"
        "Theta_deg,"
        "U1,"
        "U2,"
        "U3,"
        "X_deformed,"
        "Y_deformed,"
        "Z_deformed,"
        "U_AXIAL_Y,"
        "U_RADIAL,"
        "U_CIRCUMFERENTIAL,"
        "LE11,"
        "LE22,"
        "LE33,"
        "LE12,"
        "LE13,"
        "LE23,"
        "LE_AXIAL_Y,"
        "LE_CIRCUMFERENTIAL,"
        "LE_RADIAL,"
        "LE_YTHETA_TENSOR,"
        "GAMMA_YTHETA_ENGINEERING\n"
    )

    with open(
            path,
            "w") as output_file:

        output_file.write(header)

        for index in range(
                len(node_labels)):

            node_label = \
                node_labels[index]

            x = points[index][0]

            y = points[index][1]

            z = points[index][2]

            u1 = displacements[index][0]

            u2 = displacements[index][1]

            u3 = displacements[index][2]

            x_deformed = x + u1

            y_deformed = y + u2

            z_deformed = z + u3

            global_le = \
                global_le_values[index]

            le11 = global_le[0]

            le22 = global_le[1]

            le33 = global_le[2]

            le12 = global_le[3]

            le13 = global_le[4]

            le23 = global_le[5]

            local = local_values[index]

            output_file.write(
                "%s,%d,%.12e,"
                "%s,%d,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e,%.12e,"
                "%.12e,%.12e\n"
                % (
                    STEP_NAME,
                    frame_id,
                    frame_value,
                    instance_name,
                    node_label,
                    x,
                    y,
                    z,
                    local["RADIUS"],
                    local["THETA_RAD"],
                    local["THETA_DEG"],
                    u1,
                    u2,
                    u3,
                    x_deformed,
                    y_deformed,
                    z_deformed,
                    local["U_AXIAL_Y"],
                    local["U_RADIAL"],
                    local[
                        "U_CIRCUMFERENTIAL"
                    ],
                    le11,
                    le22,
                    le33,
                    le12,
                    le13,
                    le23,
                    local["LE_AXIAL_Y"],
                    local[
                        "LE_CIRCUMFERENTIAL"
                    ],
                    local["LE_RADIAL"],
                    local[
                        "LE_YTHETA_TENSOR"
                    ],
                    local[
                        "GAMMA_YTHETA_ENGINEERING"
                    ]
                )
            )


# ============================================================
# VTU WRITER
# ============================================================

def write_vtu(
        path,
        points,
        cells,
        point_data):
    """
    Write an ASCII VTU unstructured-grid file.
    """

    connectivity = []

    offsets = []

    cell_types = []

    offset = 0

    for cell in cells:

        connectivity.extend(cell)

        offset += len(cell)

        offsets.append(offset)

        if len(cell) == 3:

            # VTK_TRIANGLE
            cell_types.append(5)

        elif len(cell) == 4:

            # VTK_QUAD
            cell_types.append(9)

        else:

            # VTK_POLYGON
            cell_types.append(7)

    with open(
            path,
            "w") as output_file:

        output_file.write(
            '<?xml version="1.0"?>\n'
        )

        output_file.write(
            '<VTKFile '
            'type="UnstructuredGrid" '
            'version="0.1" '
            'byte_order="LittleEndian">\n'
        )

        output_file.write(
            "  <UnstructuredGrid>\n"
        )

        output_file.write(
            '    <Piece '
            'NumberOfPoints="%d" '
            'NumberOfCells="%d">\n'
            % (
                len(points),
                len(cells)
            )
        )

        # ====================================================
        # POINT DATA
        # ====================================================

        output_file.write(
            '      <PointData '
            'Scalars="LE_AXIAL_Y" '
            'Vectors="Displacement">\n'
        )

        for name, values in \
                point_data.items():

            if name == "Displacement":

                output_file.write(
                    '        <DataArray '
                    'type="Float64" '
                    'Name="%s" '
                    'NumberOfComponents="3" '
                    'format="ascii">\n'
                    % name
                )

                for value in values:

                    output_file.write(
                        "          "
                        "%.12e %.12e %.12e\n"
                        % (
                            value[0],
                            value[1],
                            value[2]
                        )
                    )

            else:

                output_file.write(
                    '        <DataArray '
                    'type="Float64" '
                    'Name="%s" '
                    'format="ascii">\n'
                    % name
                )

                for value in values:

                    output_file.write(
                        "          %.12e\n"
                        % value
                    )

            output_file.write(
                "        </DataArray>\n"
            )

        output_file.write(
            "      </PointData>\n"
        )

        # ====================================================
        # POINT COORDINATES
        # ====================================================

        output_file.write(
            "      <Points>\n"
        )

        output_file.write(
            '        <DataArray '
            'type="Float64" '
            'NumberOfComponents="3" '
            'format="ascii">\n'
        )

        for point in points:

            output_file.write(
                "          "
                "%.12e %.12e %.12e\n"
                % (
                    point[0],
                    point[1],
                    point[2]
                )
            )

        output_file.write(
            "        </DataArray>\n"
        )

        output_file.write(
            "      </Points>\n"
        )

        # ====================================================
        # CELLS
        # ====================================================

        output_file.write(
            "      <Cells>\n"
        )

        # ----------------------------------------------------
        # Connectivity
        # ----------------------------------------------------

        output_file.write(
            '        <DataArray '
            'type="Int32" '
            'Name="connectivity" '
            'format="ascii">\n'
        )

        for value in connectivity:

            output_file.write(
                "          %d\n"
                % value
            )

        output_file.write(
            "        </DataArray>\n"
        )

        # ----------------------------------------------------
        # Offsets
        # ----------------------------------------------------

        output_file.write(
            '        <DataArray '
            'type="Int32" '
            'Name="offsets" '
            'format="ascii">\n'
        )

        for value in offsets:

            output_file.write(
                "          %d\n"
                % value
            )

        output_file.write(
            "        </DataArray>\n"
        )

        # ----------------------------------------------------
        # Cell types
        # ----------------------------------------------------

        output_file.write(
            '        <DataArray '
            'type="UInt8" '
            'Name="types" '
            'format="ascii">\n'
        )

        for value in cell_types:

            output_file.write(
                "          %d\n"
                % value
            )

        output_file.write(
            "        </DataArray>\n"
        )

        output_file.write(
            "      </Cells>\n"
        )

        output_file.write(
            "    </Piece>\n"
        )

        output_file.write(
            "  </UnstructuredGrid>\n"
        )

        output_file.write(
            "</VTKFile>\n"
        )


# ============================================================
# PVD WRITER
# ============================================================

def write_pvd(
        path,
        datasets):
    """
    Write the ParaView PVD time-series file.

    datasets contains:

        (frame_time, relative_vtu_filename)
    """

    with open(
            path,
            "w") as output_file:

        output_file.write(
            '<?xml version="1.0"?>\n'
        )

        output_file.write(
            '<VTKFile '
            'type="Collection" '
            'version="0.1" '
            'byte_order="LittleEndian">\n'
        )

        output_file.write(
            "  <Collection>\n"
        )

        for time_value, filename in \
                datasets:

            output_file.write(
                '    <DataSet '
                'timestep="%.12e" '
                'group="" '
                'part="0" '
                'file="%s"/>\n'
                % (
                    time_value,
                    filename
                )
            )

        output_file.write(
            "  </Collection>\n"
        )

        output_file.write(
            "</VTKFile>\n"
        )


# ============================================================
# MAIN
# ============================================================

def main():
    """
    Main export procedure.
    """

    if not os.path.isfile(
            ODB_PATH):

        raise RuntimeError(
            "ODB file not found: %s"
            % ODB_PATH
        )

    if not os.path.exists(
            OUTPUT_FOLDER):

        os.makedirs(
            OUTPUT_FOLDER
        )

    odb = openOdb(
        path=ODB_PATH,
        readOnly=True
    )

    try:

        assembly = odb.rootAssembly

        surface, instance_name = \
            find_surface_and_instance(
                assembly
            )

        instance = \
            assembly.instances[
                instance_name
            ]

        (
            node_labels,
            points,
            cells,
            surface_element_labels
        ) = build_surface_mesh(
            surface,
            instance
        )

        frames_to_export = \
            get_frames_to_export(
                odb
            )

        print(
            "============================================================"
        )

        print(
            "ROUND-BAR OUTER CYLINDRICAL SURFACE EXPORT"
        )

        print(
            "============================================================"
        )

        print(
            "ODB:",
            ODB_PATH
        )

        print(
            "Step:",
            STEP_NAME
        )

        print(
            "Instance:",
            instance_name
        )

        print(
            "Surface:",
            SURFACE_NAME
        )

        print(
            "Surface elements:",
            len(surface_element_labels)
        )

        print(
            "Surface nodes:",
            len(node_labels)
        )

        print(
            "Surface cells:",
            len(cells)
        )

        print(
            "Frames to export:",
            len(frames_to_export)
        )

        print(
            "Output folder:",
            OUTPUT_FOLDER
        )

        print(
            "============================================================"
        )

        pvd_datasets = []

        # ====================================================
        # FRAME LOOP
        # ====================================================

        for frame_id, frame in \
                frames_to_export:

            print(
                "Processing frame %d; "
                "time = %.12e"
                % (
                    frame_id,
                    frame.frameValue
                )
            )

            # ------------------------------------------------
            # Extract and average LE
            # ------------------------------------------------

            nodal_le = \
                average_le_at_surface_nodes(
                    frame,
                    surface,
                    instance_name
                )

            # ------------------------------------------------
            # Extract displacement
            # ------------------------------------------------

            displacement_dictionary = \
                get_displacement_dictionary(
                    frame,
                    instance_name
                )

            displacements = []

            global_le_values = []

            local_values = []

            # ------------------------------------------------
            # VTU point data
            # ------------------------------------------------

            point_data = OrderedDict()

            point_data[
                "NodeID"
            ] = []

            point_data[
                "Displacement"
            ] = []

            point_data[
                "U_AXIAL_Y"
            ] = []

            point_data[
                "U_RADIAL"
            ] = []

            point_data[
                "U_CIRCUMFERENTIAL"
            ] = []

            point_data[
                "LE11"
            ] = []

            point_data[
                "LE22"
            ] = []

            point_data[
                "LE33"
            ] = []

            point_data[
                "LE12"
            ] = []

            point_data[
                "LE13"
            ] = []

            point_data[
                "LE23"
            ] = []

            point_data[
                "LE_AXIAL_Y"
            ] = []

            point_data[
                "LE_CIRCUMFERENTIAL"
            ] = []

            point_data[
                "LE_RADIAL"
            ] = []

            point_data[
                "LE_YTHETA_TENSOR"
            ] = []

            point_data[
                "GAMMA_YTHETA_ENGINEERING"
            ] = []

            point_data[
                "Radius"
            ] = []

            point_data[
                "Theta_deg"
            ] = []

            # =================================================
            # SURFACE NODE LOOP
            # =================================================

            for index, node_label in \
                    enumerate(node_labels):

                if node_label not in \
                        nodal_le:

                    raise RuntimeError(
                        "No LE value for node %d "
                        "in frame %d."
                        % (
                            node_label,
                            frame_id
                        )
                    )

                displacement = \
                    displacement_dictionary.get(
                        node_label,
                        (
                            0.0,
                            0.0,
                            0.0
                        )
                    )

                global_le = \
                    nodal_le[node_label]

                local = cylindrical_values(
                    points[index],
                    displacement,
                    global_le
                )

                displacements.append(
                    displacement
                )

                global_le_values.append(
                    global_le
                )

                local_values.append(
                    local
                )

                # --------------------------------------------
                # Node information
                # --------------------------------------------

                point_data[
                    "NodeID"
                ].append(
                    float(node_label)
                )

                # --------------------------------------------
                # Global displacement
                # --------------------------------------------

                point_data[
                    "Displacement"
                ].append(
                    displacement
                )

                # --------------------------------------------
                # Local displacement
                # --------------------------------------------
    

                # --------------------------------------------
                # Local displacement
                # --------------------------------------------

                point_data[
                    "U_AXIAL_Y"
                ].append(
                    local[
                        "U_AXIAL_Y"
                    ]
                )

                point_data[
                    "U_RADIAL"
                ].append(
                    local[
                        "U_RADIAL"
                    ]
                )

                point_data[
                    "U_CIRCUMFERENTIAL"
                ].append(
                    local[
                        "U_CIRCUMFERENTIAL"
                    ]
                )

                # --------------------------------------------
                # Global strain
                # --------------------------------------------

                point_data[
                    "LE11"
                ].append(
                    global_le[0]
                )

                point_data[
                    "LE22"
                ].append(
                    global_le[1]
                )

                point_data[
                    "LE33"
                ].append(
                    global_le[2]
                )

                point_data[
                    "LE12"
                ].append(
                    global_le[3]
                )

                point_data[
                    "LE13"
                ].append(
                    global_le[4]
                )

                point_data[
                    "LE23"
                ].append(
                    global_le[5]
                )

                # --------------------------------------------
                # Cylindrical strain
                # --------------------------------------------

                point_data[
                    "LE_AXIAL_Y"
                ].append(
                    local[
                        "LE_AXIAL_Y"
                    ]
                )

                point_data[
                    "LE_CIRCUMFERENTIAL"
                ].append(
                    local[
                        "LE_CIRCUMFERENTIAL"
                    ]
                )

                point_data[
                    "LE_RADIAL"
                ].append(
                    local[
                        "LE_RADIAL"
                    ]
                )

                point_data[
                    "LE_YTHETA_TENSOR"
                ].append(
                    local[
                        "LE_YTHETA_TENSOR"
                    ]
                )

                point_data[
                    "GAMMA_YTHETA_ENGINEERING"
                ].append(
                    local[
                        "GAMMA_YTHETA_ENGINEERING"
                    ]
                )

                # --------------------------------------------
                # Cylindrical position
                # --------------------------------------------

                point_data[
                    "Radius"
                ].append(
                    local[
                        "RADIUS"
                    ]
                )

                point_data[
                    "Theta_deg"
                ].append(
                    local[
                        "THETA_DEG"
                    ]
                )

            # =================================================
            # OUTPUT FILENAMES
            # =================================================

            vtu_filename = (
                "RB_L0_25_Y_OUTER_CYL_SURF_"
                "frame_%04d.vtu"
                % frame_id
            )

            csv_filename = (
                "RB_L0_25_Y_OUTER_CYL_SURF_"
                "frame_%04d.csv"
                % frame_id
            )

            vtu_path = os.path.join(
                OUTPUT_FOLDER,
                vtu_filename
            )

            csv_path = os.path.join(
                OUTPUT_FOLDER,
                csv_filename
            )

            # =================================================
            # WRITE VTU
            # =================================================

            write_vtu(
                vtu_path,
                points,
                cells,
                point_data
            )

            # =================================================
            # WRITE CSV
            # =================================================

            if EXPORT_CSV:

                write_csv(
                    csv_path,
                    frame_id,
                    frame.frameValue,
                    instance_name,
                    node_labels,
                    points,
                    displacements,
                    global_le_values,
                    local_values
                )

            # =================================================
            # ADD FRAME TO PVD
            # =================================================

            pvd_datasets.append(
                (
                    frame.frameValue,
                    vtu_filename
                )
            )

            print(
                "Saved VTU:",
                vtu_path
            )

            if EXPORT_CSV:

                print(
                    "Saved CSV:",
                    csv_path
                )

        # ====================================================
        # WRITE PVD TIME SERIES
        # ====================================================

        write_pvd(
            PVD_PATH,
            pvd_datasets
        )

        print(
            "============================================================"
        )

        print(
            "EXPORT FINISHED"
        )

        print(
            "PVD file:",
            PVD_PATH
        )

        print(
            "Open this PVD file in ParaView."
        )

        print(
            "============================================================"
        )

    finally:

        odb.close()


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    main()