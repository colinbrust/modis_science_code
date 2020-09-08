import subprocess as sp
import os
import ee

ee.Initialize()


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


def export_folder_to_drive(f_dir: str, out_dir: str, ee_loc: str = 'earthengine', scale: int = 500) -> None:
    """
    Exports folder of ee assets or a single asset to a google drive folder
    :param f_dir: Directory where assets are located
    :param out_dir: Google drive folder to export assets to
    :param ee_loc: string of path location of earthengine executable
    :param scale: Spatial scale to export image at
    :return: None, exports images to google drive
    """
    assets = sp.check_output([ee_loc, '--no-use_cloud_api', 'ls', f_dir]).decode().split('\n')

    for asset in assets:

        out = ee.Image(asset)

        task = ee.batch.Export.image.toDrive(
            image=out,
            description=os.path.basename(asset),
            folder=out_dir,
            scale=scale,
            maxPixels=1e13,
            region=out.geometry().bounds().getInfo()['coordinates'],
            fileNamePrefix=os.path.basename(asset)
        )

        print('Started task {}'.format(os.path.basename(asset)))
        task.start()
