# -*- coding: utf-8 -*-
"""
Check the Abaqus/Explicit ODB for the round-bar L0 = 25 mm model.

Model orientation:
    Specimen/loading axis = global Y

Important axial quantities:
    Displacement = U2
    Velocity     = V2
    Acceleration = A2
    Reaction     = RF2
    Stress       = S22
    Strain       = LE22

Run using:
    abaqus python check_RB_L0_25_Y_MBW_odb.py

The ODB and this script must be located in the same folder.
"""

from __future__ import print_function

from odbAccess import openOdb
import os
import sys


# ============================================================
# USER SETTINGS
# ============================================================

odb_name = "RB_L0_25_Y_MBW.odb"

# Explicit step duration used in the model
expected_final_time = 0.005

# Numerical tolerance for the final-time check
time_tolerance = 1.0e-6


# ============================================================
# BUILD ODB PATH
# ============================================================

script_folder = os.path.dirname(os.path.abspath(__file__))
odb_path = os.path.join(script_folder, odb_name)

if not os.path.isfile(odb_path):

    print("============================================================")
    print("ERROR: ODB FILE NOT FOUND")
    print("============================================================")
    print("Expected ODB:")
    print(odb_path)
    print("")
    print("Check that the filename is exactly:")
    print(odb_name)
    print("============================================================")

    sys.exit(1)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def get_scalar_value(field_value):
    """
    Return a scalar value from an Abaqus field-output value.

    Abaqus scalar data may be stored as:
        float
        one-component array
        one-component tuple
    """

    data = field_value.data

    try:
        return float(data)

    except (TypeError, ValueError):
        pass

    try:

        if len(data) > 0:
            return float(data[0])

    except (TypeError, ValueError, IndexError):
        pass

    return None


def get_component_value(field_value, component_index):
    """
    Extract one component from an Abaqus vector or tensor value.

    Vector ordering:
        index 0 = component 1
        index 1 = component 2
        index 2 = component 3

    Symmetric 3D tensor ordering:
        index 0 = 11
        index 1 = 22
        index 2 = 33
        index 3 = 12
        index 4 = 13
        index 5 = 23
    """

    data = field_value.data

    try:

        if len(data) > component_index:
            return float(data[component_index])

    except (TypeError, ValueError, IndexError):
        pass

    return None


def print_scalar_field_range(frame, variable_name):
    """
    Print minimum and maximum values of a scalar field.
    """

    if variable_name not in frame.fieldOutputs:

        print(variable_name + ": NOT available")
        return None

    field = frame.fieldOutputs[variable_name]
    scalar_values = []

    for field_value in field.values:

        scalar_value = get_scalar_value(field_value)

        if scalar_value is not None:
            scalar_values.append(scalar_value)

    if len(scalar_values) == 0:

        print(
            variable_name +
            ": available, but no scalar values could be extracted"
        )

        return None

    minimum_value = min(scalar_values)
    maximum_value = max(scalar_values)

    print(
        variable_name +
        ": values = " + str(len(scalar_values)) +
        ", minimum = " + str(minimum_value) +
        ", maximum = " + str(maximum_value)
    )

    return minimum_value, maximum_value


def print_component_range(
        frame,
        variable_name,
        component_index,
        component_name):
    """
    Print minimum and maximum values of one vector or tensor component.
    """

    if variable_name not in frame.fieldOutputs:

        print(component_name + ": " + variable_name + " NOT available")
        return None

    field = frame.fieldOutputs[variable_name]
    component_values = []

    for field_value in field.values:

        component_value = get_component_value(
            field_value,
            component_index
        )

        if component_value is not None:
            component_values.append(component_value)

    if len(component_values) == 0:

        print(
            component_name +
            ": could not be extracted from " +
            variable_name
        )

        return None

    minimum_value = min(component_values)
    maximum_value = max(component_values)
    maximum_absolute_value = max(
        abs(minimum_value),
        abs(maximum_value)
    )

    print(
        component_name +
        ": values = " + str(len(component_values)) +
        ", minimum = " + str(minimum_value) +
        ", maximum = " + str(maximum_value) +
        ", maximum absolute = " + str(maximum_absolute_value)
    )

    return minimum_value, maximum_value


def print_region_names(odb):
    """
    Print sets and surfaces stored at assembly and instance levels.
    """

    assembly = odb.rootAssembly

    print("")
    print("ASSEMBLY REGIONS")
    print("----------------")

    print(
        "Assembly node sets:",
        list(assembly.nodeSets.keys())
    )

    print(
        "Assembly element sets:",
        list(assembly.elementSets.keys())
    )

    print(
        "Assembly surfaces:",
        list(assembly.surfaces.keys())
    )

    print("")
    print("INSTANCE REGIONS")
    print("----------------")

    for instance_name, instance in assembly.instances.items():

        print("")
        print("Instance:", instance_name)

        try:
            print(
                "Node sets:",
                list(instance.nodeSets.keys())
            )

        except Exception:
            print("Node sets: could not be read")

        try:
            print(
                "Element sets:",
                list(instance.elementSets.keys())
            )

        except Exception:
            print("Element sets: could not be read")

        try:
            print(
                "Surfaces:",
                list(instance.surfaces.keys())
            )

        except Exception:
            print("Surfaces: could not be read")


def find_surface_by_keywords(odb, keywords):
    """
    Search assembly and instance surfaces using name keywords.

    Returns:
        region, displayed_name

    Returns:
        None, None
    when no matching surface is found.
    """

    assembly = odb.rootAssembly

    # Search assembly-level surfaces
    for surface_name, surface in assembly.surfaces.items():

        upper_name = surface_name.upper()

        for keyword in keywords:

            if keyword.upper() in upper_name:
                return surface, surface_name

    # Search instance-level surfaces
    for instance_name, instance in assembly.instances.items():

        try:
            surface_items = instance.surfaces.items()

        except Exception:
            continue

        for surface_name, surface in surface_items:

            upper_name = surface_name.upper()

            for keyword in keywords:

                if keyword.upper() in upper_name:

                    displayed_name = (
                        instance_name +
                        "." +
                        surface_name
                    )

                    return surface, displayed_name

    return None, None


def sum_field_component_on_region(
        frame,
        variable_name,
        component_index,
        region,
        displayed_region_name):
    """
    Sum one field component over a surface or node-set region.

    This is useful for summing RF2 on the loaded end surface.
    """

    if variable_name not in frame.fieldOutputs:

        print(variable_name + ": NOT available")
        return None

    try:

        subset = frame.fieldOutputs[variable_name].getSubset(
            region=region
        )

    except Exception as error:

        print(
            "Could not create " +
            variable_name +
            " subset for region " +
            displayed_region_name
        )

        print("Reason:", str(error))
        return None

    total_value = 0.0
    number_of_values = 0

    for field_value in subset.values:

        component_value = get_component_value(
            field_value,
            component_index
        )

        if component_value is not None:

            total_value += component_value
            number_of_values += 1

    if number_of_values == 0:

        print(
            "No " +
            variable_name +
            " values were found on region " +
            displayed_region_name
        )

        return None

    print(
        "Sum of " +
        variable_name +
        " component 2 on " +
        displayed_region_name +
        " = " +
        str(total_value) +
        " using " +
        str(number_of_values) +
        " values"
    )

    return total_value


def print_history_output_information(step):
    """
    Print all available history-output regions and selected nodal histories.
    """

    print("")
    print("HISTORY-OUTPUT REGIONS")
    print("----------------------")

    if len(step.historyRegions) == 0:

        print("WARNING: No history-output regions are available.")
        return

    for region_name, history_region in step.historyRegions.items():

        print("")
        print("Region:", region_name)

        available_histories = list(
            history_region.historyOutputs.keys()
        )

        print(
            "Available histories:",
            available_histories
        )

        for output_name in ["U2", "V2", "A2", "RF2"]:

            if output_name in history_region.historyOutputs:

                history_data = (
                    history_region
                    .historyOutputs[output_name]
                    .data
                )

                if len(history_data) > 0:

                    final_time = history_data[-1][0]
                    final_value = history_data[-1][1]

                    print(
                        output_name +
                        ": final time = " +
                        str(final_time) +
                        ", final value = " +
                        str(final_value)
                    )


def print_energy_histories(step):
    """
    Print final whole-model energy values.
    """

    print("")
    print("WHOLE-MODEL ENERGY HISTORIES")
    print("----------------------------")

    energy_names = [
        "ALLKE",
        "ALLIE",
        "ALLAE",
        "ALLCD",
        "ALLDMD",
        "ALLFD",
        "ALLPD",
        "ALLSE",
        "ALLVD",
        "ALLWK",
        "ALLMW",
        "ETOTAL"
    ]

    final_energy_values = {}
    energy_found = False

    for region_name, history_region in step.historyRegions.items():

        available_histories = history_region.historyOutputs

        for energy_name in energy_names:

            if energy_name in available_histories:

                energy_data = (
                    available_histories[energy_name]
                    .data
                )

                if len(energy_data) > 0:

                    energy_found = True
                    final_value = energy_data[-1][1]

                    final_energy_values[energy_name] = final_value

                    print(
                        energy_name +
                        ": final value = " +
                        str(final_value)
                    )

    if not energy_found:

        print("No whole-model energy histories were found.")
        return

    print("")
    print("ENERGY RATIOS")
    print("-------------")

    if (
            "ALLKE" in final_energy_values and
            "ALLWK" in final_energy_values):

        allke = final_energy_values["ALLKE"]
        allwk = final_energy_values["ALLWK"]

        if abs(allwk) > 0.0:

            ratio_ke_work = abs(allke / allwk)

            print(
                "|ALLKE / ALLWK| = " +
                str(ratio_ke_work)
            )

        else:

            print(
                "|ALLKE / ALLWK| could not be calculated because ALLWK = 0"
            )

    if (
            "ALLAE" in final_energy_values and
            "ALLWK" in final_energy_values):

        allae = final_energy_values["ALLAE"]
        allwk = final_energy_values["ALLWK"]

        if abs(allwk) > 0.0:

            ratio_ae_work = abs(allae / allwk)

            print(
                "|ALLAE / ALLWK| = " +
                str(ratio_ae_work)
            )

        else:

            print(
                "|ALLAE / ALLWK| could not be calculated because ALLWK = 0"
            )

    print("")
    print("ENERGY NOTE")
    print("-----------")
    print(
        "When a VUMAT does not update enerInternNew and enerInelasNew, "
        "ALLIE and ETOTAL may not represent the material energy correctly."
    )


# ============================================================
# OPEN AND CHECK ODB
# ============================================================

odb = None

try:

    odb = openOdb(
        path=odb_path,
        readOnly=True
    )

    print("============================================================")
    print("ROUND-BAR L0 = 25 mm GLOBAL-Y ODB CHECK")
    print("============================================================")
    print("ODB:", odb_path)
    print("ODB steps:", list(odb.steps.keys()))
    print(
        "Assembly instances:",
        list(odb.rootAssembly.instances.keys())
    )

    print_region_names(odb)

    if len(odb.steps) == 0:

        print("")
        print("ERROR: The ODB contains no analysis steps.")
        sys.exit(1)

    for step_name, step in odb.steps.items():

        number_of_frames = len(step.frames)

        print("")
        print("============================================================")
        print("STEP INFORMATION")
        print("============================================================")
        print("Step:", step_name)
        print("Procedure:", step.procedure)
        print("Number of frames:", number_of_frames)

        if number_of_frames == 0:

            print("WARNING: This step contains no result frames.")
            continue

        first_frame = step.frames[0]
        last_frame = step.frames[-1]

        first_time = first_frame.frameValue
        last_time = last_frame.frameValue

        print("First frame time:", first_time)
        print("Last frame time:", last_time)

        time_difference = abs(
            last_time - expected_final_time
        )

        if time_difference <= time_tolerance:

            print(
                "FINAL TIME CHECK: OK. Step reached " +
                str(expected_final_time)
            )

        else:

            print("FINAL TIME CHECK: WARNING")

            print(
                "Expected final time:",
                expected_final_time
            )

            print(
                "Actual final time:",
                last_time
            )

            print(
                "Difference:",
                time_difference
            )

        print("")
        print("LAST-FRAME FIELD OUTPUTS")
        print("------------------------")

        print(
            list(last_frame.fieldOutputs.keys())
        )

        # ====================================================
        # GLOBAL-Y AXIAL RESULTS
        # ====================================================

        print("")
        print("GLOBAL-Y AXIAL QUANTITIES")
        print("-------------------------")

        u2_range = print_component_range(
            last_frame,
            variable_name="U",
            component_index=1,
            component_name="U2"
        )

        print_component_range(
            last_frame,
            variable_name="V",
            component_index=1,
            component_name="V2"
        )

        print_component_range(
            last_frame,
            variable_name="A",
            component_index=1,
            component_name="A2"
        )

        print_component_range(
            last_frame,
            variable_name="RF",
            component_index=1,
            component_name="RF2"
        )

        print_component_range(
            last_frame,
            variable_name="S",
            component_index=1,
            component_name="S22"
        )

        print_component_range(
            last_frame,
            variable_name="LE",
            component_index=1,
            component_name="LE22"
        )

        # ====================================================
        # NOMINAL STRAIN CHECK
        # ====================================================

        print("")
        print("NOMINAL AXIAL STRAIN CHECK")
        print("--------------------------")

        if u2_range is not None:

            minimum_u2 = u2_range[0]
            maximum_u2 = u2_range[1]

            maximum_absolute_u2 = max(
                abs(minimum_u2),
                abs(maximum_u2)
            )

            nominal_strain = maximum_absolute_u2 / 25.0

            print(
                "Maximum absolute U2 = " +
                str(maximum_absolute_u2) +
                " mm"
            )

            print(
                "Nominal strain = max(|U2|) / 25 = " +
                str(nominal_strain)
            )

        # ====================================================
        # REACTION-FORCE SUM ON TOP SURFACE
        # ====================================================

        print("")
        print("TOP-END REACTION-FORCE CHECK")
        print("----------------------------")

        top_surface, top_surface_name = find_surface_by_keywords(
            odb,
            [
                "END_25V",
                "END_25",
                "TOP"
            ]
        )

        if top_surface is not None:

            sum_field_component_on_region(
                last_frame,
                variable_name="RF",
                component_index=1,
                region=top_surface,
                displayed_region_name=top_surface_name
            )

        else:

            print(
                "No top surface named END_25V, END_25 or TOP "
                "was found in the ODB."
            )

            print(
                "RF2 is still available as field output, but the script "
                "cannot automatically sum it over the top face."
            )

        # ====================================================
        # MBW STATE VARIABLES
        # ====================================================

        print("")
        print("MBW VUMAT STATE VARIABLES")
        print("-------------------------")

        important_sdvs = [
            "SDV1",    # Equivalent plastic strain
            "SDV2",    # Damage
            "SDV17",
            "SDV28",   # Active/deletion flag
            "SDV29",   # Stress triaxiality
            "SDV30",   # Lode parameter
            "SDV31",   # Damage-initiation indicator
            "SDV32",   # Ductile failure indicator
            "SDV33",
            "SDV34",
            "SDV36",   # Cleavage failure indicator
            "SDV37",
            "SDV38",
            "SDV39",
            "SDV40",
            "SDV41",
            "SDV42"
        ]

        for sdv_name in important_sdvs:

            print_scalar_field_range(
                last_frame,
                sdv_name
            )

        # ====================================================
        # ELEMENT STATUS
        # ====================================================

        print("")
        print("ELEMENT STATUS")
        print("--------------")

        status_range = print_scalar_field_range(
            last_frame,
            "STATUS"
        )

        if status_range is not None:

            status_minimum = status_range[0]
            status_maximum = status_range[1]

            if (
                    status_minimum == 1.0 and
                    status_maximum == 1.0):

                print(
                    "STATUS CHECK: All elements are active."
                )

            else:

                print(
                    "STATUS CHECK: Some elements may be deleted or inactive."
                )

        # ====================================================
        # GENERAL OUTPUT AVAILABILITY
        # ====================================================

        print("")
        print("GENERAL FIELD-OUTPUT AVAILABILITY")
        print("---------------------------------")

        general_variables = [
            "S",
            "LE",
            "U",
            "V",
            "A",
            "RF",
            "COORD",
            "STATUS"
        ]

        for variable_name in general_variables:

            if variable_name in last_frame.fieldOutputs:

                number_of_values = len(
                    last_frame
                    .fieldOutputs[variable_name]
                    .values
                )

                print(
                    variable_name +
                    ": available, values = " +
                    str(number_of_values)
                )

            else:

                print(
                    variable_name +
                    ": NOT available"
                )

        # ====================================================
        # HISTORY OUTPUTS AND ENERGIES
        # ====================================================

        print_history_output_information(step)
        print_energy_histories(step)

        print("")
        print("============================================================")
        print("END OF STEP CHECK")
        print("============================================================")

finally:

    if odb is not None:
        odb.close()


print("")
print("============================================================")
print("RB_L0_25_Y_MBW ODB CHECK FINISHED")
print("============================================================")