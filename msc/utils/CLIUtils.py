import subprocess as sp


def stop_all_tasks(ee_loc: str = 'earthengine') -> None:
    """
    :param ee_loc: string of path location of earthengine executable.
    :return: None, stops all tasks that are running in Earthengine
    """
    task_list = sp.check_output([ee_loc, 'task', 'list']).decode('utf-8').split('\n')
    for task in task_list:

        components = ' '.join(task.split()).split()

        if len(components) > 0:
            task_id = components[0]

            if 'READY' in components or 'RUNNING' in components:
                sp.call([ee_loc, '--no-use_cloud_api', 'task', 'cancel', task_id])


def task_running(ee_loc: str = 'earthengine') -> bool:
    """
    :param ee_loc: string of path location of earthengine executable.
    :return: Boolean stating whether or not any tasks are running in earthengine.
    """
    out = sp.check_output([ee_loc, 'task', 'list']).decode()

    if (out.find('READY') != -1) or (out.find('RUNNING') != -1):
        return True
    else:
        return False
