let sessionActive = true;

async function addMessage(message, isUser = false) {
    const messagesDiv = document.getElementById('chat-messages');
    const messageDiv = document.createElement('div');
    messageDiv.className = `message ${isUser ? 'user-message' : 'bot-message'}`;

    // Create message header
    const headerDiv = document.createElement('div');
    headerDiv.className = 'message-header';
    
    const avatarDiv = document.createElement('div');
    avatarDiv.className = 'message-avatar';
    avatarDiv.style.background = isUser ? 'var(--message-user-bg)' : 'var(--accent-color)';
    avatarDiv.textContent = isUser ? 'You' : 'AI';
    
    const nameSpan = document.createElement('span');
    nameSpan.textContent = isUser ? 'You' : 'Assistant';
    
    headerDiv.appendChild(avatarDiv);
    headerDiv.appendChild(nameSpan);
    messageDiv.appendChild(headerDiv);
    
    // Create message content
    const contentDiv = document.createElement('div');
    contentDiv.className = 'message-content';
    
    // Format the message based on the response structure
    if (isUser) {
        contentDiv.textContent = message;
    } else {
        // Format the response data
        let formattedMessage = `Query: ${message.query}\n`;
        formattedMessage += `Time period: ${message.timeframe}\n`;
        formattedMessage += `Total papers found: ${message.total_papers}\n`;
        formattedMessage += `Total trials found: ${message.total_trials}\n\n`;
        
        // Add metrics summary
        formattedMessage += "Papers Metrics:\n";
        for (const [category, items] of Object.entries(message.metrics_summary.papers)) {
            formattedMessage += `${category}: ${items.join(', ')}\n`;
        }
        
        formattedMessage += "\nTrials Metrics:\n";
        for (const [category, items] of Object.entries(message.metrics_summary.trials)) {
            formattedMessage += `${category}: ${items.join(', ')}\n`;
        }
        
        formattedMessage += "\nAnalysis:\n";
        formattedMessage += message.analysis;
        
        contentDiv.innerHTML = marked.parse(formattedMessage);  // Use marked to render markdown
    }
    
    messageDiv.appendChild(contentDiv);
    messagesDiv.appendChild(messageDiv);
    
    // Smooth scroll to bottom
    messagesDiv.scrollTo({
        top: messagesDiv.scrollHeight,
        behavior: 'smooth'
    });
}

function showTypingIndicator() {
    document.getElementById('typing-indicator').style.display = 'block';
    const messagesDiv = document.getElementById('chat-messages');
    messagesDiv.scrollTo({
        top: messagesDiv.scrollHeight,
        behavior: 'smooth'
    });
}

function hideTypingIndicator() {
    document.getElementById('typing-indicator').style.display = 'none';
}

async function sendQuery() {
    const queryInput = document.getElementById('query-input');
    const sendButton = document.querySelector('.send-button');
    const query = queryInput.value.trim();
    
    if (!query) return;

    // Disable input while processing
    queryInput.disabled = true;
    sendButton.disabled = true;

    addMessage(query, true);
    queryInput.value = '';
    
    showTypingIndicator();

    try {
        const response = await fetch('/process-query', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ query: query })
        });

        const data = await response.json();
        hideTypingIndicator();
        
        if (response.ok) {
            addMessage(data, false);  // Pass the entire data object
        } else {
            addMessage(`Error: ${data.detail}`, false);
        }
    } catch (error) {
        hideTypingIndicator();
        addMessage('Error processing query: ' + error.message, false);
    } finally {
        // Re-enable input
        queryInput.disabled = false;
        sendButton.disabled = false;
        queryInput.focus();
    }
}

async function generateReport() {
    const queryInput = document.getElementById('query-input');
    const generateButton = document.querySelector('.generate-button');
    const query = queryInput.value.trim();
    
    if (!query) return;

    // Show loading state
    document.getElementById('loading').style.display = 'block';
    document.getElementById('report-section').style.display = 'none';
    generateButton.disabled = true;

    try {
        const response = await fetch('/process-query', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ query: query })
        });

        const data = await response.json();
        console.log("Received data:", data); // Debug log
        
        if (response.ok) {
            displayReport(data);
        } else {
            alert(`Error: ${data.detail}`);
        }
    } catch (error) {
        console.error("Error:", error); // Debug log
        alert('Error generating report: ' + error.message);
    } finally {
        document.getElementById('loading').style.display = 'none';
        generateButton.disabled = false;
    }
}

function displayReport(data) {
    // Show report section
    const reportSection = document.getElementById('report-section');
    reportSection.style.display = 'block';

    // Set report title
    document.getElementById('report-title').textContent = 
        `Analysis Report: ${data.query}`;

    // Display overview metrics
    document.getElementById('overview-metrics').innerHTML = `
        <p><strong>Time Period:</strong> ${data.timeframe}</p>
        <p><strong>Total Papers:</strong> ${data.total_papers}</p>
        <p><strong>Total Trials:</strong> ${data.total_trials}</p>
    `;

    // Display paper metrics
    document.getElementById('paper-metrics').innerHTML = formatMetrics(data.metrics_summary.papers);

    // Display trial metrics
    document.getElementById('trial-metrics').innerHTML = formatMetrics(data.metrics_summary.trials);

    // Display detailed analysis
    document.getElementById('detailed-analysis').innerHTML = marked.parse(data.analysis);
}

function formatMetrics(metrics) {
    return Object.entries(metrics)
        .map(([category, items]) => `
            <div class="metric-category">
                <h4>${category.charAt(0).toUpperCase() + category.slice(1)}</h4>
                <ul>
                    ${items.map(item => `<li>${item}</li>`).join('')}
                </ul>
            </div>
        `).join('');
}

// Handle Enter key in input
document.getElementById('query-input').addEventListener('keypress', function(e) {
    if (e.key === 'Enter' && !e.shiftKey) {
        e.preventDefault();
        generateReport();
    }
});

// Clean up when window closes
window.addEventListener('beforeunload', async () => {
    if (sessionActive) {
        try {
            await fetch('/cleanup', { method: 'POST' });
        } catch (error) {
            console.error('Error during cleanup:', error);
        }
    }
});