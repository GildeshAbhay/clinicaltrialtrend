import os
from workflows.logs.log import setup_logger

logger = setup_logger()

def list_pdf_files(pdf_directory):
    pdf_files_list = []  # List to store tuples of (filename, full_path)
    logger.info(f"Listing PDF files in {pdf_directory}")
    if not os.path.exists(pdf_directory):
        print(f"Directory '{pdf_directory}' does not exist!")
        return pdf_files_list
    
    try:
        # Walk through directory and all subdirectories
        for root, dirs, files in os.walk(pdf_directory):
            for file in files:
                if file.lower().endswith('.pdf'):
                    full_path = os.path.join(root, file)
                    pdf_files_list.append( full_path)
        logger.info(f"Found {len(pdf_files_list)} PDF files")
        return pdf_files_list
                
    except Exception as e:
        print(f"Error accessing directory: {e}")
        logger.error(f"Error accessing directory: {e}")
        return pdf_files_list
