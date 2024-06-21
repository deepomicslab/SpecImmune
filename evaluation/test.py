import requests
from bs4 import BeautifulSoup

# URL of the PharmVar download page
url = 'https://www.pharmvar.org/download'

# Step 1: Send a GET request to the download page
response = requests.get(url)
response.raise_for_status()  # Check that the request was successful

# Step 2: Parse the HTML content using BeautifulSoup
soup = BeautifulSoup(response.text, 'html.parser')

# Step 3: Find the download link (assuming it's in an <a> tag with a specific class)
download_link = soup.find('a', class_='action-button-download-complete')

if download_link and 'href' in download_link.attrs:
    file_url = download_link['href']
    if not file_url.startswith('http'):
        file_url = 'https://www.pharmvar.org' + file_url

    # Step 4: Download the file
    file_response = requests.get(file_url)
    file_response.raise_for_status()  # Check that the request was successful

    # Step 5: Save the file locally
    with open('complete_database.zip', 'wb') as file:
        file.write(file_response.content)
    
    print('Download complete.')
else:
    print('Download link not found.')