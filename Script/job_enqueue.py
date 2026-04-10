import sys
import json
import redis
import subprocess
from rq import Queue, Worker

# Connect to Redis
redis_conn = redis.Redis()
queue = Queue('bio_jobs', connection=redis_conn, default_timeout=10800)

# Get job data from PHP
try:
    job_data = json.loads(sys.argv[1])
except json.JSONDecodeError as e:
    print(f"JSON Decode Error: {e}")
    sys.exit(1)

# Enqueue the job
job = queue.enqueue("Script.worker.process_job", job_data)
print(f"Job Enqueued: ID={job.id}, Status={job.get_status()}")

# Check if at least one worker is running
workers = Worker.all(connection=redis_conn)

if len(workers) == 0:  # Start worker only if no worker is running
    print("No active workers found. Starting a worker...")
    subprocess.Popen(["rq", "worker", "bio_jobs"],
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
