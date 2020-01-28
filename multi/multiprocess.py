import concurrent.futures

from log_setup import setup_logger

logger = setup_logger(__name__, 'logs/multi.log')

MAX_PROCESS_NUMBER = 8
MAX_JOB_NUMBER = MAX_PROCESS_NUMBER


def launch_multiprocess(runs, func):
    """
    Spawn multiprocesses avoiding memory problems.
    The 'runs' parameter must be a list of 'cmd' s.t., subprocess.Popen(cmd).wait() works.
    :param runs: a list of cmd (parameters, executable, etc).
    :param fn: the function to be parallelized
    :return:
    """
    done = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_PROCESS_NUMBER) as executor:
        jobs = {}
        runs_left = len(runs)
        runs_iter = iter(runs)

        while runs_left:
            logger.info('Multiprocess {}: runs_left = {}'.format(func, runs_left))
            for run in runs_iter:
                job = executor.submit(func, *run)
                jobs[job] = run
                if len(jobs) > MAX_JOB_NUMBER:
                   break

            for job in concurrent.futures.as_completed(jobs):
                result = job.result()
                runs_left -= 1
                run = jobs[job]
                del jobs[job]
                done.append(result)
                break

    logger.info('Multiprocess {} ended'.format(func))
    return done

