import concurrent.futures

from log_setup import setup_logger

logger = setup_logger(__name__, 'logs/multi.log')

MAX_THREAD_NUMBER = 8
MAX_JOB_NUMBER = 8


def launch_multithreads(runs, func):
    logger.info('Starting multithreaded {} '.format(func))

    done = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_THREAD_NUMBER) as executor:
        jobs = {}
        runs_left = len(runs)
        runs_iter = iter(runs)

        while runs_left:
            logger.info('Multithread {}: runs_left = {}'.format(func, runs_left))
            for run in runs_iter:
                job = executor.submit(func, *run)
                jobs[job] = run
                if len(jobs) > MAX_JOB_NUMBER:
                    break

            for job in concurrent.futures.as_completed(jobs):
                runs_left -= 1
                res = job.result()
                run = jobs[job]
                logger.debug('Ended job {}'.format(run))
                del jobs[job]
                done.append(res)
                break

    logger.info('Multithreaded {} ended'.format(func))
    return done