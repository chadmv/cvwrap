import maya.cmds as cmds
import maya.mel as mel
import maya.OpenMayaUI as OpenMayaUI
import os

NAME_WIDGET = 'cvwrap_name'
RADIUS_WIDGET = 'cvwrap_radius'
NEW_BIND_MESH_WIDGET = 'cvwrap_newbindmesh'
BIND_FILE_WIDGET = 'cvwrap_bindfile'
MENU_ITEMS = []

def create_menuitems():
    global MENU_ITEMS
    if MENU_ITEMS:
        # Already created
        return
    for menu in ['mainDeformMenu', 'mainRigDeformationsMenu']:
        # Make sure the menu widgets exist first.
        mel.eval('ChaDeformationsMenu MayaWindow|{0};'.format(menu))
        items = cmds.menu(menu, q=True, ia=True)
        for item in items:
            if cmds.menuItem(item, q=True, divider=True):
                section = cmds.menuItem(item, q=True, label=True)
            if cmds.menuItem(item, q=True, label=True) == 'Wrap':
                if section == 'Create':
                    cvwrap_item = cmds.menuItem(label="cvWrap", command=create_cvwrap,
                                                sourceType='python', insertAfter=item, parent=menu)
                    cvwrap_options = cmds.menuItem(command=display_cvwrap_options,
                                                   insertAfter=cvwrap_item, parent=menu,
                                                   optionBox=True)
                    MENU_ITEMS.append(cvwrap_item)
                    MENU_ITEMS.append(cvwrap_options)
                elif section == 'Edit':
                    cvwrap_item = cmds.menuItem(label="cvWrap Binding",
                                                command=display_binding_dialog,
                                                sourceType='python', insertAfter=item, parent=menu)
                    MENU_ITEMS.append(cvwrap_item)


def create_cvwrap(*args, **kwargs):
    cmds.loadPlugin('cvwrap', qt=True)
    sel = cmds.ls(sl=True)
    if len(sel) >= 2:
        kwargs = get_command_kwargs()
        result = cmds.cvWrap(**kwargs)
        print result
    else:
        raise RuntimeError("Select at least one surface and one influence object.")

def get_command_kwargs():
    """Gets the cvWrap command arguments either from the option box widgets or the saved
    option vars.  If the widgets exist, their values will be saved to the option vars.
    @return A dictionary of the kwargs to the cvWrap command."""
    args = {}
    if cmds.textFieldGrp(NAME_WIDGET, exists=True):
        args['name'] = cmds.textFieldGrp(NAME_WIDGET, q=True, text=True)
        cmds.optionVar(sv=(NAME_WIDGET, args['name']))
    else:
        args['name'] = cmds.optionVar(q=NAME_WIDGET) or 'cvWrap#'
    if cmds.floatSliderGrp(RADIUS_WIDGET, exists=True):
        args['radius'] = cmds.floatSliderGrp(RADIUS_WIDGET, q=True, value=True)
        cmds.optionVar(fv=(RADIUS_WIDGET, args['radius']))
    else:
        args['radius'] = cmds.optionVar(q=RADIUS_WIDGET)

    if cmds.checkBoxGrp(NEW_BIND_MESH_WIDGET, exists=True):
        if cmds.checkBoxGrp(NEW_BIND_MESH_WIDGET, q=True, v1=True):
            args['newBindMesh'] = True
            cmds.optionVar(iv=(NEW_BIND_MESH_WIDGET, 1))
        else:
            cmds.optionVar(iv=(NEW_BIND_MESH_WIDGET, 0))
    else:
        value = cmds.optionVar(q=NEW_BIND_MESH_WIDGET)
        if value:
            args['newBindMesh'] = True

    if cmds.textFieldButtonGrp(BIND_FILE_WIDGET, exists=True):
        bind_file = cmds.textFieldButtonGrp(BIND_FILE_WIDGET, q=True, text=True)
        bind_file = os.path.expandvars(bind_file.strip())
        if bind_file:
            if os.path.exists(bind_file):
                args['binding'] = bind_file
            else:
                cmds.warning('{0} does not exist.'.format(bind_file))

    return args


def display_cvwrap_options(*args, **kwargs):
    cmds.loadPlugin('cvwrap', qt=True)
    layout = mel.eval('getOptionBox')
    cmds.setParent(layout)
    cmds.columnLayout(adj=True)

    for widget in [NAME_WIDGET, RADIUS_WIDGET, BIND_FILE_WIDGET, NEW_BIND_MESH_WIDGET]:
        # Delete the widgets so we don't create multiple controls with the same name
        try:
            cmds.deleteUI(widget, control=True)
        except:
            pass

    cmds.textFieldGrp(NAME_WIDGET, label='Node name', text='cvWrap#')
    radius = cmds.optionVar(q=RADIUS_WIDGET)
    cmds.floatSliderGrp(RADIUS_WIDGET, label='Sample radius', field=True, minValue=0.0,
                        maxValue=100.0, fieldMinValue=0.0, fieldMaxValue=100.0, value=radius,
                        step=0.01, precision=2)
    cmds.textFieldButtonGrp(BIND_FILE_WIDGET, label='Binding file ', text='', buttonLabel='Browse',
                            bc=display_bind_file_dialog)
    use_new_bind_mesh = cmds.optionVar(q=NEW_BIND_MESH_WIDGET)
    cmds.checkBoxGrp(NEW_BIND_MESH_WIDGET, numberOfCheckBoxes=1, label='Create new bind mesh',
                     v1=use_new_bind_mesh)
    mel.eval('setOptionBoxTitle("cvWrap Options");')
    mel.eval('setOptionBoxCommandName("cvWrap");')
    apply_close_button = mel.eval('getOptionBoxApplyAndCloseBtn;')
    cmds.button(apply_close_button, e=True, command=apply_and_close)
    apply_button = mel.eval('getOptionBoxApplyBtn;')
    cmds.button(apply_button, e=True, command=create_cvwrap)
    reset_button = mel.eval('getOptionBoxResetBtn;')
    # For some reason, the buttons in the menu only accept MEL.
    cmds.button(reset_button, e=True,
                command='python("import cvwrap.menu; cvwrap.menu.reset_to_defaults()");')
    close_button = mel.eval('getOptionBoxCloseBtn;')
    cmds.button(close_button, e=True, command=close_option_box)
    save_button = mel.eval('getOptionBoxSaveBtn;')
    cmds.button(save_button, e=True,
                command='python("import cvwrap.menu; cvwrap.menu.get_command_kwargs()");')
    mel.eval('showOptionBox')


def apply_and_close(*args, **kwargs):
    """Create the cvWrap deformer and close the option box."""
    create_cvwrap()
    mel.eval('saveOptionBoxSize')
    close_option_box()


def close_option_box(*args, **kwargs):
    mel.eval('hideOptionBox')


def display_bind_file_dialog(*args, **kwargs):
    """Displays the dialog to choose the binding file with which to create the cvWrap deformer."""
    root_dir = cmds.workspace(q=True, rootDirectory=True)
    start_directory = os.path.join(root_dir, 'data')
    file_path = cmds.fileDialog2(fileFilter='*.wrap', dialogStyle=2, fileMode=1,
                                 startingDirectory=start_directory)
    if file_path:
        cmds.textFieldButtonGrp(BIND_FILE_WIDGET, e=True, text=file_path[0])


def reset_to_defaults(*args, **kwargs):
    """Reset the cvWrap option box widgets to their defaults."""
    cmds.textFieldGrp(NAME_WIDGET, e=True, text='cvWrap#')
    cmds.floatSliderGrp(RADIUS_WIDGET, e=True, value=0)
    cmds.textFieldButtonGrp(BIND_FILE_WIDGET, e=True, text='')
    cmds.checkBoxGrp(NEW_BIND_MESH_WIDGET, e=True, v1=False)


def display_binding_dialog(*args, **kwargs):
    cmds.loadPlugin('cvwrap', qt=True)


def destroy_menuitems():
    """Remove the cvWrap items from the menus."""
    global MENU_ITEMS
    for item in MENU_ITEMS:
        cmds.deleteUI(item, menuItem=True)
    MENU_ITEMS = []

