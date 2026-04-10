import os
import subprocess
import mysql.connector
import datetime
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders
import zipfile

def zip_directory(source_dir, zip_name, exclude_file="PSSM_Features_ML_17W.csv"):
    """
    Compress all files and folders inside `source_dir` into a zip file,
    excluding the specified file.
    
    :param source_dir: Directory to zip.
    :param zip_name: Name of the output zip file (without extension).
    :param exclude_file: File to exclude from the zip.
    """
    zip_path = f"{zip_name}.zip"
    
    # Create the directory for the zip file if it doesn't exist.
    os.makedirs(os.path.dirname(zip_path), exist_ok=True)
    
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                if file == exclude_file:
                    continue  # Skip the excluded file
                file_path = os.path.join(root, file)
                archive_name = os.path.relpath(file_path, source_dir)
                zipf.write(file_path, archive_name)

    print(f"✅ Successfully created zip: {zip_path}")
# Email Configuration
SMTP_SERVER = "smtp.gmail.com"  # Change if using another email provider
SMTP_PORT = 587
EMAIL_SENDER = "mayankk3549@gmail.com"
EMAIL_PASSWORD = "sjcorbqfjkntlvyk"
SENDER_NAME="S8kPred"
def send_email1(job_id, user_email, status):
    """Send an email notification about job completion."""
    subject = f"Job {job_id} Status: {'Completed' if status == 2 else 'Failed'}"
    body = f"Hello,\n\nYour job (ID: {job_id}) has been {'successfully completed' if status == 2 else 'failed'}.\n\nBest Regards,\nYour Team"
    
    msg = MIMEText(body)
    msg["Subject"] = subject
    msg["From"] = EMAIL_SENDER
    msg["To"] = user_email

    try:
        server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
        server.starttls()
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, user_email, msg.as_string())
        server.quit()
        print(f"Email sent to {user_email} for Job {job_id}.")
    except Exception as e:
        print(f"Failed to send email: {str(e)}")
        
        
def send_email(job_id, user_email, status, attachment_path=None):
    """Send an email notification with an HTML template and an optional attachment."""
    
    # Define email subject & message
    subject = f"Job {job_id} Status: {'Completed' if status == 2 else 'Failed'}"
    
    html_content = f"""
    <html>
    <body style="font-family: Arial, sans-serif; color: #333;">
        <!-- <h2 style="color: #007BFF;">Job Notification</h2> -->
        <p>Hello,</p>
        <p>Your job <b>ID: {job_id}</b> has been <b style="color:{'green' if status == 2 else 'red'};">
        {'successfully completed' if status == 2 else 'failed'}</b>.</p>
        <!-- <p>Please check the attached results (if available).</p> -->
       <p>Click the button below to view your job results:</p>
		<a href='https://s8kpred.in//result.php?JobID={job_id}' class='btn'>View Job Results</a>
        <p class='footer'>If you have any questions, please contact support at <a href='mailto:mayank2801@gmail.com'>mayank2801@gmail.com</a>.</p>
        <hr>
        <p>Best Regards,<br><b>S8kPred Team</b></p>
    </body>
    </html>
    """

    # Create email message
    msg = MIMEMultipart()
    msg["From"] = f"{SENDER_NAME} <{EMAIL_SENDER}>"
    msg["To"] = user_email
    msg["Subject"] = subject

    # Attach HTML message
    msg.attach(MIMEText(html_content, "html"))

    # Attach file if provided and exists
    if attachment_path and os.path.exists(attachment_path):
        with open(attachment_path, "rb") as attachment:
            part = MIMEBase("application", "octet-stream")
            part.set_payload(attachment.read())
            encoders.encode_base64(part)
            part.add_header("Content-Disposition", f"attachment; filename={os.path.basename(attachment_path)}")
            msg.attach(part)
    else:
        print(f"Attachment not found at: {attachment_path}")

    # Send email
    try:
        server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
        server.starttls()
        server.login(EMAIL_SENDER, EMAIL_PASSWORD)
        server.sendmail(EMAIL_SENDER, user_email, msg.as_string())
        server.quit()
        print(f"✅ Email sent to {user_email} for Job {job_id}.")
    except Exception as e:
        print(f"❌ Failed to send email: {str(e)}")
        
def update_job_status_run(job_id, status, user_email):
    """Update job status in MySQL and send an email notification."""
    try:
        conn = mysql.connector.connect(
            host="localhost",
            user="s8kpred",
            password="S8kpred@1234",
            database="s8kpred"
        )
        cursor = conn.cursor()

        # Update job status and completion time
        query = "UPDATE jobstatus SET JobStatus=%s WHERE JobID=%s"
        cursor.execute(query, (status, job_id))

        conn.commit()
        cursor.close()
        conn.close()

        print(f"Job {job_id} status updated to {status}.")

    except Exception as e:
        print(f"Error updating job status: {str(e)}")
def update_job_status(job_id, status, user_email):
    """Update job status in MySQL and send an email notification."""
    try:
        conn = mysql.connector.connect(
            host="localhost",
            user="s8kpred",
            password="S8kpred@1234",
            database="s8kpred"
        )
        cursor = conn.cursor()

        # Update job status and completion time
        query = "UPDATE jobstatus SET JobStatus=%s, JobCompletionTime=%s WHERE JobID=%s"
        completion_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        cursor.execute(query, (status, completion_time, job_id))

        conn.commit()
        cursor.close()
        conn.close()

        print(f"Job {job_id} status updated to {status}.")
        #zip_directory("Jobs/{job_id}", "Jobs/{job_id}/result")
        #attachment_path="Jobs/{job_id}/result.zip" 
        # Send email to user
        send_email(job_id, user_email, status)

    except Exception as e:
        print(f"Error updating job status: {str(e)}")

def process_job(job_data):
    job_id = job_data['job_id']
    fasta_path = job_data['fasta_path']
    option = job_data['secondary_structure_option']
    user_email = job_data['user_email']  # Get user email from job data

    job_dir = f"Jobs/{job_id}"
    log_file = os.path.join(job_dir, "log.dat")
    #log=open(log_file,mode='a')
    pssm_output_dir = os.path.join(job_dir, "pssm_outputs")

    try:
        # **STEP 1: Generate PSSM**
        update_job_status_run(job_id, 1, user_email)
        pssm_cmd = f"python Script/PSSMBasedFeaturesGeneration.py generate_pssm_files {job_id}"
        with open(log_file, "a") as log:
            log.write(f"Running: {pssm_cmd}\n")
        subprocess.run(pssm_cmd, shell=True, check=True)

        # **Check if PSSM file was created**
        if not any(fname.endswith(".pssm") for fname in os.listdir(pssm_output_dir)):
            with open(log_file, "a") as log:
                log.write("No PSSM File Generated\n")
            update_job_status(job_id, 3, user_email)  # FAILED
            return

        # **STEP 2: Run secondary structure prediction**
        if option == "ThreeSecondaryStructureOption":
            ss_cmd = f"python Script/SecondaryStructurePredictionThreeState.py ThreeStateSSPred {job_id}"
            expected_output = os.path.join(job_dir, "ResultThreeState.fas")

        elif option == "EightSecondaryStructureOption":
            ss_cmd = f"python Script/SecondaryStructurePredictionEightState.py EightStateSSPred {job_id}"
            expected_output = os.path.join(job_dir, "ResultEightState.fas")

        elif option == "BothSecondaryStructureOption":
            ss_cmd1 = f"python Script/SecondaryStructurePredictionThreeState.py ThreeStateSSPred {job_id}"
            ss_cmd2 = f"python Script/SecondaryStructurePredictionEightState.py EightStateSSPred {job_id}"
            expected_output1 = os.path.join(job_dir, "ResultThreeState.fas")
            expected_output2 = os.path.join(job_dir, "ResultEightState.fas")
            with open(log_file, "a") as log:
                log.write(f"Running: {ss_cmd1}\n")
            subprocess.run(ss_cmd1, shell=True, check=True)
            if not os.path.exists(expected_output1):
                with open(log_file, "a") as log:
                    log.write("Unable to predict Three State Secondary Structure\n")
                update_job_status(job_id, 3, user_email)  # FAILED
                return
            with open(log_file, "a") as log:
                log.write(f"Running: {ss_cmd2}\n")
            subprocess.run(ss_cmd2, shell=True, check=True)
            if not os.path.exists(expected_output2):
                with open(log_file, "a") as log:
                    log.write("Unable to predict Eight State Secondary Structure\n")
                update_job_status(job_id, 3, user_email)  # FAILED
                return
            with open(log_file, "a") as log:
                log.write(f"Sucessfully Job Completed\n")
            update_job_status(job_id, 2, user_email)  # SUCCESS
            return
        with open(log_file, "a") as log:
            log.write(f"Running: {ss_cmd}\n")
        subprocess.run(ss_cmd, shell=True, check=True)

        # **Check if output file was created**
        if not os.path.exists(expected_output):
            with open(log_file, "a") as log:
                log.write("Unable to predict Secondary Structure\n")
            update_job_status(job_id, 3, user_email)  # FAILED
            return

        update_job_status(job_id, 2, user_email)  # SUCCESS

    except subprocess.CalledProcessError:
        update_job_status(job_id, 3, user_email)  # FAILED
